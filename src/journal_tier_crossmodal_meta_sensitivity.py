"""
Cross-modal primary claim × two-cohort DL meta (exploratory sensitivity).

When outputs/exploratory_meta_DE_two_cohort_DL.csv exists (GSE188646 + GSE87102),
recomputes Spearman(miRNA weighted target burden, aging logFC) using:
  (i) GSE188646 pseudobulk logFC (cohort1; same as journal_tier_crossmodal), and
  (ii) DerSimonian–Laird pooled beta_DL on the same symbol axis.

Does not assert exchangeability across sex/assay — documents as convergent-context check only.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from config import OUTPUT_DIR


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    r = stats.spearmanr(x, y, nan_policy="omit")
    stat = float(getattr(r, "statistic", getattr(r, "correlation", np.nan)))
    pv = float(getattr(r, "pvalue", np.nan))
    return stat, pv


def run_crossmodal_meta_sensitivity(out_dir: Path | None = None, log=print) -> None:
    out_dir = out_dir or OUTPUT_DIR
    burden_p = out_dir / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv"
    meta_p = out_dir / "exploratory_meta_DE_two_cohort_DL.csv"
    if not burden_p.is_file():
        log("Cross-modal meta sensitivity: missing burden CSV; skipped.")
        return
    if not meta_p.is_file():
        log(
            "Cross-modal meta sensitivity: missing exploratory_meta_DE_two_cohort_DL.csv "
            "(cohort2 limma or DL meta failed, or cohort1 missing se_logFC — see SA_COMPLETENESS_CHECK.txt); skipped."
        )
        return

    b = pd.read_csv(burden_p)
    m = pd.read_csv(meta_p)
    if "gene" not in b.columns or "weighted_burden" not in b.columns:
        log("Cross-modal meta sensitivity: unexpected burden columns; skipped.")
        return
    if "gene" not in m.columns or "beta_DL" not in m.columns:
        log("Cross-modal meta sensitivity: unexpected meta columns; skipped.")
        return

    b["gene"] = b["gene"].astype(str).str.strip().str.upper()
    m["gene"] = m["gene"].astype(str).str.strip().str.upper()
    log_col = "logFC" if "logFC" in b.columns else None
    if log_col is None:
        log("Cross-modal meta sensitivity: burden CSV missing logFC; skipped.")
        return

    j = b.merge(m, on="gene", how="inner", suffixes=("_burden", "_meta"))
    j = j.dropna(subset=["weighted_burden", log_col, "beta_DL"])
    if len(j) < 200:
        log(f"Cross-modal meta sensitivity: too few genes after merge ({len(j)}); skipped.")
        return

    r_c1, p_c1 = _spearman(j["weighted_burden"].values.astype(float), j[log_col].values.astype(float))
    r_dl, p_dl = _spearman(j["weighted_burden"].values.astype(float), j["beta_DL"].values.astype(float))
    r_c2 = float("nan")
    p_c2 = float("nan")
    if "logFC_gse87102" in j.columns:
        r_c2, p_c2 = _spearman(j["weighted_burden"].values.astype(float), j["logFC_gse87102"].values.astype(float))

    out_csv = out_dir / "exploratory_crossmodal_burden_vs_two_cohort_meta.csv"
    keep = ["gene", "weighted_burden", "n_mirnas", log_col, "beta_DL", "se_DL"]
    if "logFC_gse188646" in j.columns:
        keep.append("logFC_gse188646")
    if "logFC_gse87102" in j.columns:
        keep.append("logFC_gse87102")
    j[[c for c in keep if c in j.columns]].to_csv(out_csv, index=False)

    summ = {
        "completed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_genes_intersection": int(len(j)),
        "spearman_rho_burden_vs_gse188646_logFC": r_c1,
        "spearman_p_burden_vs_gse188646_logFC": p_c1,
        "spearman_rho_burden_vs_meta_beta_DL": r_dl,
        "spearman_p_burden_vs_meta_beta_DL": p_dl,
        "spearman_rho_burden_vs_gse87102_logFC": r_c2,
        "spearman_p_burden_vs_gse87102_logFC": p_c2,
        "methodology_note": "Sensitivity for primary cross-modal coupling: same miRNA burden vector vs cohort1 logFC vs two-study DL pooled logFC. Sex/assay differ between cohorts; meta beta is contextual, not a causal estimand.",
        "gene_universe_note": "Rows = intersection of crossmodal burden CSV with genes in two-cohort meta CSV (n differs from exploratory_crossmodal_mirna_aging_summary.json n_genes_merged).",
        "output_genes_csv": str(out_csv.resolve()),
    }
    (out_dir / "exploratory_crossmodal_meta_cohort_sensitivity_summary.json").write_text(
        json.dumps(summ, indent=2), encoding="utf-8"
    )
    log(f"Cross-modal meta sensitivity: wrote {out_csv.name} and exploratory_crossmodal_meta_cohort_sensitivity_summary.json")
    log(json.dumps({k: summ[k] for k in summ if k != "methodology_note"}, indent=2))
