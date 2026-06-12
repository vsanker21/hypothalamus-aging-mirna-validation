"""
Allen MERFISH (ABC-WMB) spatial anatomy validation for GSE188646 strata.

Uses MERFISH cell metadata with Allen CCF parcellation (no expression matrix download).
Maps GSE188646 strata via HypoMap Spearman reference (if available) to MERFISH WMB subclasses
(Tanycyte NN, Ependymal NN, etc.) and quantifies hypothalamus / V3-adjacent CCF coordinates.

Outputs:
  exploratory_merfish_spatial_validation_per_stratum.csv
  exploratory_merfish_tanycyte_ccf_coordinates_summary.csv
  exploratory_merfish_spatial_validation_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from abc_atlas_fetch import download_abc_file

MERFISH_NICHE_SUBCLASS_MAP = {
    "Tanycytes": "Tanycyte NN",
    "Tanycyte": "Tanycyte NN",
    "Ependymal": "Ependymal NN",
    "ParsTuber": "Tanycyte NN",
    "Astroependymal": "Astroependymal NN",
    "Hypendymal": "Hypendymal NN",
}
HY_PARCELLATION_DIVISION = "HY"
V3_CCF_BOUNDS = {
    "x_ccf_min": 5.5,
    "x_ccf_max": 8.5,
    "y_ccf_min": 6.0,
    "y_ccf_max": 9.5,
    "z_ccf_min": 3.0,
    "z_ccf_max": 7.0,
}


def _load_merfish_parcellation(log) -> pd.DataFrame | None:
    path = download_abc_file("merfish_cell_metadata_parcellation", log=log)
    if path is None or not path.is_file():
        return None
    usecols = [
        "cell_label",
        "subclass",
        "supertype",
        "parcellation_division",
        "parcellation_structure",
        "parcellation_substructure",
        "x_ccf",
        "y_ccf",
        "z_ccf",
    ]
    try:
        return pd.read_csv(path, usecols=usecols, low_memory=False)
    except Exception as exc:
        if log:
            log(f"MERFISH parcellation read failed: {exc}")
        return None


def _in_v3_bounds(df: pd.DataFrame) -> pd.Series:
    return (
        (df["x_ccf"].astype(float) >= V3_CCF_BOUNDS["x_ccf_min"])
        & (df["x_ccf"].astype(float) <= V3_CCF_BOUNDS["x_ccf_max"])
        & (df["y_ccf"].astype(float) >= V3_CCF_BOUNDS["y_ccf_min"])
        & (df["y_ccf"].astype(float) <= V3_CCF_BOUNDS["y_ccf_max"])
        & (df["z_ccf"].astype(float) >= V3_CCF_BOUNDS["z_ccf_min"])
        & (df["z_ccf"].astype(float) <= V3_CCF_BOUNDS["z_ccf_max"])
    )


def _hypomap_to_merfish_subclass(ref_type: str) -> str | None:
    rt = str(ref_type).lower()
    if "ependymal" in rt:
        return "ependymal"
    if "tanycyte" in rt or "parstuber" in rt:
        return "tanycyte"
    if "hypendymal" in rt:
        return "hypendymal"
    if "astroependymal" in rt:
        return "astroependymal"
    return None


def run_merfish_spatial_validation(
    out_dir: Path,
    log,
    *,
    rho_threshold: float = 0.82,
) -> None:
    log("\n=== Allen MERFISH spatial anatomy validation (ABC-WMB CCF parcellation) ===")
    meta = _load_merfish_parcellation(log)
    if meta is None or meta.empty:
        log("MERFISH validation: parcellation metadata unavailable; skipped.")
        return

    niche_subs = sorted(set(MERFISH_NICHE_SUBCLASS_MAP.values()))
    niche_re = r"Tanycyte|Ependymal|Hypendymal|Astroependymal"
    niche_cells = meta[meta["subclass"].astype(str).str.contains(niche_re, case=False, na=False, regex=True)].copy()
    hy_cells = meta[meta["parcellation_division"].astype(str) == HY_PARCELLATION_DIVISION].copy()

    ccf_rows = []
    for sub in sorted(niche_cells["subclass"].astype(str).unique()):
        sub_cells = niche_cells[niche_cells["subclass"].astype(str) == sub]
        sub_hy = sub_cells[sub_cells["parcellation_division"] == HY_PARCELLATION_DIVISION]
        v3_mask = _in_v3_bounds(sub_cells)
        ccf_rows.append(
            {
                "merfish_subclass": sub,
                "n_cells_total": len(sub_cells),
                "n_cells_hypothalamus_division": len(sub_hy),
                "frac_hypothalamus": len(sub_hy) / len(sub_cells) if len(sub_cells) else float("nan"),
                "n_cells_v3_ccf_bounds": int(v3_mask.sum()),
                "frac_v3_ccf_bounds": float(v3_mask.mean()) if len(sub_cells) else float("nan"),
                "median_x_ccf": float(sub_cells["x_ccf"].median()) if len(sub_cells) else float("nan"),
                "median_y_ccf": float(sub_cells["y_ccf"].median()) if len(sub_cells) else float("nan"),
                "median_z_ccf": float(sub_cells["z_ccf"].median()) if len(sub_cells) else float("nan"),
                "median_x_ccf_hypothalamus_only": float(sub_hy["x_ccf"].median()) if len(sub_hy) else float("nan"),
            }
        )
    ccf_df = pd.DataFrame(ccf_rows)
    ccf_path = out_dir / "exploratory_merfish_tanycyte_ccf_coordinates_summary.csv"
    ccf_df.to_csv(ccf_path, index=False)

    # Per-stratum: link HypoMap best type -> MERFISH subclass -> spatial stats
    hm_path = out_dir / "gse188646_hypomap_mapping" / "hypomap_custom_ref_spearman.csv"
    niche_path = out_dir / "exploratory_niche_lability_per_stratum.csv"
    rows = []
    if hm_path.is_file():
        hm = pd.read_csv(hm_path)
        hm["seurat_cluster_id"] = hm["seurat_cluster_id"].astype(str).str.strip()
        for sid, g in hm.groupby("seurat_cluster_id"):
            best = g.sort_values("rho", ascending=False).iloc[0]
            mer_sub = _hypomap_to_merfish_subclass(str(best["ref_celltype"]))
            spatial = {}
            if mer_sub:
                sc = ccf_df[ccf_df["merfish_subclass"].astype(str).str.contains(mer_sub, case=False, na=False)]
                if not sc.empty:
                    # aggregate across matching subclass rows
                    spatial = {
                        "n_cells_total": int(sc["n_cells_total"].sum()),
                        "frac_hypothalamus": float(
                            sc["n_cells_hypothalamus_division"].sum() / max(sc["n_cells_total"].sum(), 1)
                        ),
                        "frac_v3_ccf_bounds": float(
                            sc["n_cells_v3_ccf_bounds"].sum() / max(sc["n_cells_total"].sum(), 1)
                        ),
                        "median_x_ccf_hypothalamus_only": float(sc["median_x_ccf_hypothalamus_only"].median()),
                    }
            rows.append(
                {
                    "stratum": sid,
                    "hypomap_rank1_type": best["ref_celltype"],
                    "hypomap_rank1_rho": float(best["rho"]),
                    "merfish_mapped_subclass": mer_sub,
                    "merfish_n_cells_total": spatial.get("n_cells_total"),
                    "merfish_frac_hypothalamus": spatial.get("frac_hypothalamus"),
                    "merfish_frac_v3_ccf_bounds": spatial.get("frac_v3_ccf_bounds"),
                    "merfish_median_x_ccf_hy": spatial.get("median_x_ccf_hypothalamus_only"),
                    "merfish_spatial_validated": bool(
                        mer_sub
                        and spatial.get("frac_hypothalamus", 0)
                        and float(spatial.get("frac_hypothalamus", 0)) > 0.05
                    ),
                }
            )
    per = pd.DataFrame(rows)
    if niche_path.is_file():
        nd = pd.read_csv(niche_path)
        nd["stratum"] = nd["stratum"].astype(str)
        per = nd.merge(per, on="stratum", how="left")

    out_csv = out_dir / "exploratory_merfish_spatial_validation_per_stratum.csv"
    per.to_csv(out_csv, index=False)

    n_spatial = int(per.get("merfish_spatial_validated", pd.Series(dtype=bool)).astype(bool).sum())
    n_conc = 0
    if "is_third_ventricle_niche" in per.columns:
        n_conc = int(
            (per["is_third_ventricle_niche"].astype(bool) & per["merfish_spatial_validated"].astype(bool)).sum()
        )

    summ = {
        "reference": "Allen ABC-WMB MERFISH C57BL6J-638850 + CCF parcellation (adult anatomy)",
        "method": (
            "HypoMap transcriptomic best-match -> MERFISH WMB subclass -> CCF coordinates "
            "in hypothalamus (HY) division and heuristic third-ventricle CCF bounds."
        ),
        "n_merfish_cells_total": int(len(meta)),
        "n_merfish_niche_subclass_cells": int(len(niche_cells)),
        "n_hypothalamus_division_cells": int(len(hy_cells)),
        "merfish_subclass_spatial_summary": ccf_rows,
        "v3_ccf_bounds_heuristic": V3_CCF_BOUNDS,
        "n_strata_merfish_spatial_validated": n_spatial,
        "n_concordant_marker_and_merfish_spatial": n_conc,
        "caveat": (
            "MERFISH is adult-only (no aging). Spatial validation uses taxonomy mapping + CCF "
            "parcellation, not per-stratum expression correlation (ABC example-gene panel too sparse)."
        ),
    }
    json_path = out_dir / "exploratory_merfish_spatial_validation_summary.json"
    json_path.write_text(json.dumps(summ, indent=2, allow_nan=False), encoding="utf-8")
    log(f"Wrote {ccf_path.name}; MERFISH niche cells n={len(niche_cells)}")
    log(f"Wrote {out_csv.name}; spatial-validated strata n={n_spatial}, concordant n={n_conc}")
    log(f"Wrote {json_path.name}")
