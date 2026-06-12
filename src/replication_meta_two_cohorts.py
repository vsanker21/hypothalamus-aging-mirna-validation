"""
Pre-specified two-study random-effects (DerSimonian–Laird) meta-analysis on logFC + SE.

Inputs:
  outputs/gse188646_young_vs_aged_deg.csv — cohort1 (edgeR pseudobulk; must include se_logFC)
  outputs/cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv — cohort2 (limma)

Output:
  outputs/exploratory_meta_DE_two_cohort_DL.csv
  outputs/exploratory_META_TWO_COHORT_DL_README.txt

See data/provenance/REPLICATION_AND_META_POLICY.txt and data/replication/META_TWO_COHORTS_README.txt.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
from scipy.stats import norm


def _dl_two_study(b1: float, se1: float, b2: float, se2: float) -> tuple[float, float, float]:
    """DerSimonian–Laird tau^2 and pooled beta, SE(beta) for k=2 studies (inverse-variance weights)."""
    if se1 <= 0 or se2 <= 0 or not math.isfinite(b1) or not math.isfinite(b2):
        return float("nan"), float("nan"), float("nan")
    w1, w2 = 1.0 / (se1**2), 1.0 / (se2**2)
    wsum = w1 + w2
    if wsum <= 0:
        return float("nan"), float("nan"), float("nan")
    b_fe = (w1 * b1 + w2 * b2) / wsum
    q = w1 * (b1 - b_fe) ** 2 + w2 * (b2 - b_fe) ** 2
    df_ = 1.0
    c = wsum - (w1**2 + w2**2) / wsum
    tau2 = max(0.0, (q - df_) / c) if c > 1e-12 else 0.0
    v = tau2 + 1.0 / wsum
    se = math.sqrt(v) if v > 0 else float("nan")
    return b_fe, se, tau2


def run_two_cohort_meta(project_root: Path, out_dir: Path, log) -> None:
    c1 = out_dir / "gse188646_young_vs_aged_deg.csv"
    c2 = out_dir / "cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv"
    if not c1.is_file() or not c2.is_file():
        log("Meta two-cohort: need cohort1 and cohort2 DE CSVs in outputs/; skipped.")
        return
    d1 = pd.read_csv(c1)
    d2 = pd.read_csv(c2)
    if "se_logFC" not in d1.columns:
        log("Meta two-cohort: cohort1 missing se_logFC; rerun pseudobulk_edgeR_gse188646.R.")
        return
    if "se_logFC" not in d2.columns:
        log("Meta two-cohort: cohort2 missing se_logFC; skipped.")
        return
    d1["gene"] = d1["gene"].astype(str).str.upper()
    d2["gene"] = d2["gene"].astype(str).str.upper()
    m = d1.merge(d2, on="gene", suffixes=("_gse188646", "_gse87102"), how="inner")
    if len(m) < 1000:
        log(f"Meta two-cohort: only {len(m)} genes in intersection; check symbol harmonisation.")
    rows = []
    for _, r in m.iterrows():
        b1, s1 = float(r["logFC_gse188646"]), float(r["se_logFC_gse188646"])
        b2, s2 = float(r["logFC_gse87102"]), float(r["se_logFC_gse87102"])
        beta, se, tau2 = _dl_two_study(b1, s1, b2, s2)
        if math.isfinite(beta) and math.isfinite(se) and se > 0:
            z = beta / se
            p = float(2 * norm.sf(abs(z)))
        else:
            z, p = float("nan"), float("nan")
        rows.append(
            {
                "gene": r["gene"],
                "beta_DL": beta,
                "se_DL": se,
                "tau2_DL": tau2,
                "z": z,
                "p_two_sided": p,
                "logFC_gse188646": b1,
                "logFC_gse87102": b2,
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        from statsmodels.stats.multitest import multipletests

        mask = out["p_two_sided"].astype(float).notna()
        q = pd.Series(index=out.index, dtype=float)
        if mask.sum():
            _, fdr, _, _ = multipletests(out.loc[mask, "p_two_sided"].astype(float), method="fdr_bh")
            q.loc[mask] = fdr
        out["fdr_bh"] = q
    out_path = out_dir / "exploratory_meta_DE_two_cohort_DL.csv"
    out.to_csv(out_path, index=False)
    readme = out_dir / "exploratory_META_TWO_COHORT_DL_README.txt"
    readme.write_text(
        "Two-cohort DerSimonian–Laird meta on logFC + SE (exploratory).\n"
        "- Cohort1: GSE188646 pseudobulk edgeR QLF (female snRNA; Aged vs Young).\n"
        "- Cohort2: GSE87102 C57BL/6 hypothalamus microarray (male bulk; old vs young titles).\n"
        "Sex and assay differ — interpret as convergent/divergent context, not replication of identical design.\n"
        "tau^2 with k=2 is unstable; use effect directions and CI width as primary readouts.\n",
        encoding="utf-8",
    )
    log(f"Meta two-cohort: wrote {out_path.name} ({len(out)} genes).")
