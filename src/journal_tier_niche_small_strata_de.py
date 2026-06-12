"""
Targeted stratified DE for small V3-adjacent GSE188646 clusters (48, 50, 57).

Clusters 48 (Tanycyte_ependymal), 50 (Radial_glia_like), 57 (Pericyte) fall below the
default GSE188646_STRATUM_MIN_CELLS=400 threshold. Runs edgeR pseudobulk with relaxed
min cells (default 50) when RDS is available.

Outputs:
  gse188646_strata/stratum_48_young_vs_aged_deg.csv (etc.)
  exploratory_gse188646_small_niche_strata_de_summary.json
  exploratory_gse188646_small_niche_strata_de.csv
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pandas as pd

from data_acquisition import ensure_gse188646_rds_download, resolve_gse188646_rds_path

SMALL_NICHE_STRATA = ("48", "50", "57")
STRATUM_LABELS = {
    "48": "Tanycyte_ependymal",
    "50": "Radial_glia_like",
    "57": "Pericyte",
}


def _run_targeted_r(rds: Path, log, *, min_cells: int = 50) -> bool:
    proj = Path(__file__).resolve().parents[1]
    r_script = proj / "r" / "pseudobulk_stratified_edgeR_gse188646.R"
    if not r_script.is_file():
        log(f"Small-niche strata DE: R script missing {r_script}")
        return False
    env = os.environ.copy()
    env["GSE188646_STRATUM_MIN_CELLS"] = str(min_cells)
    env["GSE188646_STRATUM_MIN_GENES"] = "50"
    env["GSE188646_STRATUM_MIN_REP"] = "2"
    env["GSE188646_STRATUM_ONLY"] = ",".join(SMALL_NICHE_STRATA)
    env["GSE188646_STRATUM_APPEND_MANIFEST"] = "1"
    cmd = ["Rscript", str(r_script), str(rds), "seurat_clusters"]
    log(
        f"Small-niche strata DE: R edgeR for clusters {','.join(SMALL_NICHE_STRATA)} "
        f"(min_cells={min_cells})..."
    )
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            cwd=str(proj),
            capture_output=True,
            text=True,
            timeout=7200,
            encoding="utf-8",
            errors="replace",
            env=env,
        )
        if proc.returncode != 0:
            log(f"Small-niche strata DE failed (exit {proc.returncode})")
            if proc.stderr:
                log(proc.stderr[-3000:])
            return False
        if proc.stdout:
            log(proc.stdout.strip()[-500:])
        return True
    except FileNotFoundError:
        log("Small-niche strata DE: Rscript not found")
        return False
    except subprocess.TimeoutExpired:
        log("Small-niche strata DE: timed out")
        return False


def _summarize_stratum_de(out_dir: Path, sid: str) -> dict:
    p = out_dir / "gse188646_strata" / f"stratum_{sid}_young_vs_aged_deg.csv"
    label = STRATUM_LABELS.get(sid, sid)
    if not p.is_file():
        return {
            "stratum": sid,
            "label": label,
            "deg_available": False,
        }
    deg = pd.read_csv(p)
    padj_col = next((c for c in ("p_val_adj", "FDR", "padj") if c in deg.columns), None)
    n_fdr05 = 0
    top_genes: list[str] = []
    if padj_col:
        sig = deg[deg[padj_col].astype(float) <= 0.05].sort_values(padj_col)
        n_fdr05 = int(len(sig))
        top_genes = sig["gene"].astype(str).head(10).tolist() if "gene" in sig.columns else []
    union_path = out_dir / "mirna_target_union_genes.csv"
    n_union_sig = 0
    if union_path.is_file() and padj_col and n_fdr05:
        u = pd.read_csv(union_path)
        ucol = "gene" if "gene" in u.columns else u.columns[0]
        union = {str(g).strip().upper() for g in u[ucol].dropna()}
        sig_genes = {str(g).upper() for g in deg.loc[deg[padj_col].astype(float) <= 0.05, "gene"]}
        n_union_sig = len(sig_genes & union)
    manifest = out_dir / "gse188646_strata" / "manifest.csv"
    n_cells = None
    if manifest.is_file():
        m = pd.read_csv(manifest)
        row = m[m["stratum"].astype(str) == sid]
        if not row.empty:
            n_cells = int(row.iloc[0]["n_cells"])
    return {
        "stratum": sid,
        "label": label,
        "deg_available": True,
        "deg_path": p.name,
        "n_cells": n_cells,
        "n_genes_tested": int(len(deg)),
        "n_fdr05": n_fdr05,
        "n_fdr05_in_mirna_union": n_union_sig,
        "top_fdr05_genes": top_genes,
    }


def run_small_niche_strata_de(out_dir: Path, log) -> None:
    log("\n=== GSE188646 small niche strata DE (clusters 48/50/57) ===")
    rds = resolve_gse188646_rds_path()
    if rds is None:
        log("  RDS not cached; attempting AUTO_FETCH_GSE188646_RDS download...")
        ensure_gse188646_rds_download(log)
        rds = resolve_gse188646_rds_path()
    ran_r = False
    if rds is not None:
        ran_r = _run_targeted_r(rds, log)
    else:
        log("  GSE188646 RDS unavailable; cannot compute small-cluster pseudobulk DE.")

    rows = [_summarize_stratum_de(out_dir, sid) for sid in SMALL_NICHE_STRATA]
    df = pd.DataFrame(rows)
    csv_path = out_dir / "exploratory_gse188646_small_niche_strata_de.csv"
    df.to_csv(csv_path, index=False)

    summ = {
        "reference": "GSE188646 snRNA-seq stratified pseudobulk edgeR (Zhang et al. 2022)",
        "target_strata": list(SMALL_NICHE_STRATA),
        "stratum_labels": STRATUM_LABELS,
        "rds_used": str(rds) if rds else None,
        "edgeR_ran": ran_r,
        "min_cells_threshold": 50,
        "rationale": (
            "Clusters 48 (Tanycyte_ependymal), 50 (Radial_glia_like), and 57 (Pericyte) are "
            "V3-adjacent niche types with n=202/158/82 cells — below default min_cells=400. "
            "Relaxed threshold enables exploratory young-vs-aged DE with full replicate design."
        ),
        "per_stratum": rows,
        "caveat": "Small cell counts inflate variance; interpret as exploratory unless FDR hits replicate.",
    }
    json_path = out_dir / "exploratory_gse188646_small_niche_strata_de_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {csv_path.name}; edgeR ran={ran_r}")
    log(f"Wrote {json_path.name}")
