"""
Cross-cohort sensitivity of |effect| Mann–Whitney tests (union vs non-union).

Documents whether GSE188646 pseudobulk distributional differences in |logFC| extend to
GSE87102 microarray and DL meta |beta|. Direction (attenuated vs amplified |effect| in union)
is reported via one-sided Mann–Whitney tails.

Outputs:
  exploratory_crosscohort_lability_replication.csv
  exploratory_crosscohort_lability_replication_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def _union_set(long_path: Path) -> set[str]:
    long_df = pd.read_csv(long_path)
    tcol = "target_gene" if "target_gene" in long_df.columns else long_df.columns[-1]
    return {str(g).strip().upper() for g in long_df[tcol].astype(str)}


def _lability_stats(effect: pd.Series, genes: pd.Series, union_set: set[str]) -> dict:
    g = genes.astype(str).str.strip().str.upper()
    eff = effect.astype(float)
    in_u = g.isin(union_set)
    u = np.abs(eff.loc[in_u].values)
    o = np.abs(eff.loc[~in_u].values)
    if len(u) < 5 or len(o) < 5:
        return {}
    med_u, med_o = float(np.median(u)), float(np.median(o))
    try:
        mw = stats.mannwhitneyu(u, o, alternative="two-sided")
        mw_p = float(mw.pvalue)
        mw_p_greater = float(stats.mannwhitneyu(u, o, alternative="greater").pvalue)
        mw_p_less = float(stats.mannwhitneyu(u, o, alternative="less").pvalue)
        n1, n2 = len(u), len(o)
        r_rb = 1.0 - (2.0 * mw.statistic) / (n1 * n2)
    except ValueError:
        mw_p = mw_p_greater = mw_p_less = float("nan")
        r_rb = float("nan")
    return {
        "n_genes_union": int(in_u.sum()),
        "n_genes_nonunion": int((~in_u).sum()),
        "median_abs_effect_union": med_u,
        "median_abs_effect_nonunion": med_o,
        "delta_median_abs_effect": med_u - med_o,
        "mannwhitney_abs_effect_union_vs_nonunion_p": mw_p,
        "mannwhitney_union_greater_p": mw_p_greater,
        "mannwhitney_union_less_p": mw_p_less,
        "rank_biserial_abs_effect": float(r_rb) if np.isfinite(r_rb) else float("nan"),
    }


def run_crosscohort_lability_replication(out_dir: Path, log) -> None:
    log("\n=== Cross-cohort sensitivity: |effect| Mann–Whitney (union vs non-union) ===")
    long_path = out_dir / "mirna_targets_long.csv"
    if not long_path.is_file():
        log("Cross-cohort lability: missing mirna_targets_long.csv; skipped.")
        return
    union_set = _union_set(long_path)

    cohorts: list[tuple[str, Path, str, str]] = [
        (
            "GSE188646_pseudobulk",
            out_dir / "gse188646_young_vs_aged_deg.csv",
            "gene",
            "logFC",
        ),
        (
            "GSE87102_microarray",
            out_dir / "cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv",
            "gene",
            "logFC",
        ),
        (
            "two_cohort_DL_meta",
            out_dir / "exploratory_meta_DE_two_cohort_DL.csv",
            "gene",
            "beta_DL",
        ),
    ]

    rows: list[dict] = []
    for label, path, gcol, ecol in cohorts:
        if not path.is_file():
            log(f"Cross-cohort lability: missing {path.name}; skipped {label}.")
            continue
        df = pd.read_csv(path)
        if gcol not in df.columns or ecol not in df.columns:
            log(f"Cross-cohort lability: unexpected columns in {path.name}; skipped {label}.")
            continue
        st = _lability_stats(df[ecol], df[gcol], union_set)
        if not st:
            log(f"Cross-cohort lability: insufficient genes for {label}; skipped.")
            continue
        rows.append(
            {
                "cohort_label": label,
                "effect_column": ecol,
                "source_csv": path.name,
                "sex_assay_note": {
                    "GSE188646_pseudobulk": "Female snRNA-seq pseudobulk; edgeR QLF Aged vs Young",
                    "GSE87102_microarray": "Male hypothalamus microarray; limma Aged vs Young",
                    "two_cohort_DL_meta": "DerSimonian–Laird pooled beta; sex/assay heterogeneous",
                }.get(label, ""),
                **st,
            }
        )

    if not rows:
        log("Cross-cohort lability: no cohorts computed; skipped.")
        return

    out_csv = out_dir / "exploratory_crosscohort_lability_replication.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)

    gse = next((r for r in rows if r["cohort_label"] == "GSE188646_pseudobulk"), None)
    c2 = next((r for r in rows if r["cohort_label"] == "GSE87102_microarray"), None)
    meta = next((r for r in rows if r["cohort_label"] == "two_cohort_DL_meta"), None)

    summ = {
        "n_cohorts_tested": len(rows),
        "cohort_labels": [r["cohort_label"] for r in rows],
        "gse188646_mannwhitney_two_sided_p": gse["mannwhitney_abs_effect_union_vs_nonunion_p"] if gse else None,
        "gse188646_mannwhitney_union_less_p": gse.get("mannwhitney_union_less_p") if gse else None,
        "gse87102_mannwhitney_two_sided_p": c2["mannwhitney_abs_effect_union_vs_nonunion_p"] if c2 else None,
        "meta_dl_mannwhitney_two_sided_p": meta["mannwhitney_abs_effect_union_vs_nonunion_p"] if meta else None,
        "gse188646_union_attenuated_abs_effect": bool(
            gse and gse.get("delta_median_abs_effect", 0) < 0
        ),
        "methodology_note": (
            "Same miRTarBase target union; Mann–Whitney on |logFC| or |beta_DL| for union vs non-union "
            "within each cohort DE table. One-sided 'less' tail tests attenuated |effect| in union targets."
        ),
        "caveat": (
            "Cohorts differ by sex, assay, and tissue processing. GSE188646 shows significant attenuation "
            "of |logFC| in union genes; GSE87102 and meta do not reproduce that distributional shift."
        ),
    }
    json_path = out_dir / "exploratory_crosscohort_lability_replication_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {out_csv.name} ({len(rows)} cohorts).")
    log(f"Wrote {json_path.name}")
    for r in rows:
        log(
            f"  {r['cohort_label']}: MW P={r['mannwhitney_abs_effect_union_vs_nonunion_p']:.4g}, "
            f"delta median |effect|={r['delta_median_abs_effect']:.4f}"
        )
