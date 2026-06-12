"""
Allen Mouse Brain Atlas ISH — supplementary marker anatomy (Gpr50, Foxj1, Rax).

Queries brain-map.org RMA API for structure-level expression energy in hypothalamus /
ventricular regions using SectionDataSet → StructureUnionize (Mouse Brain product id=1).

Outputs:
  exploratory_allen_ish_marker_anatomy.csv
  exploratory_allen_ish_marker_anatomy_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ISH_MARKERS = ["Gpr50", "Foxj1", "Rax", "Tm4sf1", "Aqp4"]
STRUCTURE_QUERIES = {
    "hypothalamus": 1097,
    "third_ventricle": 135,
    "median_eminence": 476,
    "ventral_hypothalamic_zone": 467,  # VHZ — V3-adjacent hypothalamus
}
MOUSE_BRAIN_PRODUCT_ID = 1
MOUSE_ORGANISM_ID = 2
RMA_BASE = "http://api.brain-map.org/api/v2/data/query.json"


def _rma_query(criteria: str) -> list[dict]:
    params = {"criteria": criteria, "num_rows": "all"}
    try:
        r = requests.get(RMA_BASE, params=params, timeout=90)
        r.raise_for_status()
        data = r.json()
        msg = data.get("msg", [])
        if isinstance(msg, str):
            return []
        return msg if isinstance(msg, list) else []
    except Exception:
        return []


def _mouse_gene_id(acronym: str) -> int | None:
    rows = _rma_query(
        f"model::Gene,rma::criteria,[acronym$eq'{acronym}'],[organism_id$eq{MOUSE_ORGANISM_ID}]"
    )
    if not rows:
        return None
    # Prefer gene record with Mouse Brain ISH SectionDataSets (handles alias duplicates).
    for row in rows:
        if not isinstance(row, dict) or not row.get("id"):
            continue
        gid = int(row["id"])
        if _section_datasets_for_gene(gid):
            return gid
    first = rows[0]
    if isinstance(first, dict) and first.get("id"):
        return int(first["id"])
    return None


def _section_datasets_for_gene(gene_id: int) -> list[dict]:
    return _rma_query(
        f"model::SectionDataSet,rma::criteria,genes[id$eq{gene_id}],products[id$eq{MOUSE_BRAIN_PRODUCT_ID}]"
    )


def _structure_energy(section_data_set_id: int, structure_id: int) -> dict | None:
    crit = (
        f"model::StructureUnionize,rma::criteria,"
        f"[section_data_set_id$eq{section_data_set_id}],[structure_id$eq{structure_id}]"
    )
    rows = _rma_query(crit)
    if not rows:
        return None
    row = rows[0]
    if not isinstance(row, dict):
        return None
    return {
        "expression_energy": float(row["expression_energy"]) if row.get("expression_energy") is not None else None,
        "expression_density": float(row["expression_density"]) if row.get("expression_density") is not None else None,
    }


def run_allen_ish_marker_anatomy(out_dir: Path, log) -> None:
    log("\n=== Allen Mouse Brain Atlas ISH marker anatomy (supplementary) ===")
    rows = []
    for gene in ISH_MARKERS:
        gid = _mouse_gene_id(gene)
        if gid is None:
            for struct_name, struct_id in STRUCTURE_QUERIES.items():
                rows.append(
                    {
                        "gene": gene,
                        "gene_id": None,
                        "structure": struct_name,
                        "structure_id": struct_id,
                        "max_expression_energy": None,
                        "max_expression_density": None,
                        "n_section_datasets": 0,
                        "note": "gene_not_found",
                    }
                )
            continue

        sds_list = _section_datasets_for_gene(gid)
        for struct_name, struct_id in STRUCTURE_QUERIES.items():
            energies: list[float] = []
            densities: list[float] = []
            best_sds = None
            for sds in sds_list:
                sid = sds.get("id")
                if sid is None:
                    continue
                stats = _structure_energy(int(sid), struct_id)
                if stats is None or stats["expression_energy"] is None:
                    continue
                energies.append(stats["expression_energy"])
                if stats["expression_density"] is not None:
                    densities.append(stats["expression_density"])
                if best_sds is None or stats["expression_energy"] >= max(energies):
                    best_sds = int(sid)
            rows.append(
                {
                    "gene": gene,
                    "gene_id": gid,
                    "structure": struct_name,
                    "structure_id": struct_id,
                    "max_expression_energy": max(energies) if energies else None,
                    "max_expression_density": max(densities) if densities else None,
                    "n_section_datasets": len(sds_list),
                    "best_section_data_set_id": best_sds,
                }
            )

    df = pd.DataFrame(rows)
    csv_path = out_dir / "exploratory_allen_ish_marker_anatomy.csv"
    df.to_csv(csv_path, index=False)

    enrich_notes = []
    for gene in ISH_MARKERS:
        sub = df[(df["gene"] == gene) & df["structure"].notna()]
        v3 = sub.loc[sub["structure"] == "third_ventricle", "max_expression_energy"]
        hyp = sub.loc[sub["structure"] == "hypothalamus", "max_expression_energy"]
        vhz = sub.loc[sub["structure"] == "ventral_hypothalamic_zone", "max_expression_energy"]
        if len(v3) and len(hyp) and pd.notna(v3.iloc[0]) and pd.notna(hyp.iloc[0]) and hyp.iloc[0] > 0:
            enrich_notes.append(
                {
                    "gene": gene,
                    "v3_to_hypothalamus_energy_ratio": float(v3.iloc[0] / hyp.iloc[0]),
                    "v3_energy": float(v3.iloc[0]),
                    "hypothalamus_energy": float(hyp.iloc[0]),
                }
            )
        if len(v3) and len(vhz) and pd.notna(v3.iloc[0]) and pd.notna(vhz.iloc[0]) and vhz.iloc[0] > 0:
            enrich_notes.append(
                {
                    "gene": gene,
                    "v3_to_vhz_energy_ratio": float(v3.iloc[0] / vhz.iloc[0]),
                }
            )

    n_with_energy = int(df["max_expression_energy"].notna().sum())
    summ = {
        "reference": "Allen Mouse Brain Atlas ISH (adult); brain-map.org RMA API",
        "api_method": "SectionDataSet(genes, product=1) → StructureUnionize per structure_id",
        "markers_queried": ISH_MARKERS,
        "structures": STRUCTURE_QUERIES,
        "n_rows": len(df),
        "n_rows_with_expression_energy": n_with_energy,
        "v3_enrichment_heuristic": enrich_notes,
        "evidentiary_weight": (
            "low — adult ISH anatomy only (not aging or single-cell). "
            "Quantitative structure unionize energies support qualitative V3/hypothalamus localization."
        ),
        "status": "ok" if n_with_energy >= len(ISH_MARKERS) else "partial",
    }
    json_path = out_dir / "exploratory_allen_ish_marker_anatomy_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {csv_path.name} ({len(df)} rows; {n_with_energy} with energy)")
    log(f"Wrote {json_path.name} (status={summ['status']})")
