"""
Cell-type / cluster–resolved cross-modal coupling (GSE188646).

Uses per-stratum pseudobulk DE tables from r/pseudobulk_stratified_edgeR_gse188646.R
(outputs/gse188646_strata/manifest.csv + stratum_*_young_vs_aged_deg.csv) when present.

For each stratum with sufficient genes after merge with the global miRNA-target burden table:
  - Spearman(weighted_burden, logFC) on the same burden construction as journal_tier_crossmodal
  - Gene-label permutation null: shuffle logFC within the merged gene set (destroys rank coupling
    while preserving marginal logFC distribution for that stratum)
  - Mann–Whitney |logFC| in target-union vs non-union within that stratum's DE table

Exploratory only: burden is still defined from the htNSC-biased MOESM × miRTarBase program;
stratum-specific logFC is pseudobulk young vs aged within that cluster.

Enable stratified DE first, e.g.:
  GSE188646_STRATUM_COL=seurat_clusters  (or GSE188646_AUTO_STRATIFIED_PSEUDOBULK=1 in run_extended)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    r, p = stats.spearmanr(x, y, nan_policy="omit")
    if not np.isfinite(r):
        return float("nan"), float("nan")
    return float(r), float(p)


def _build_gene_burden(long_df: pd.DataFrame, wmap: dict[str, float]) -> pd.DataFrame:
    long_df = long_df.copy()
    long_df["mirna_w"] = long_df["mirna"].map(wmap).fillna(0.0).clip(lower=0.0)
    g_count = long_df.groupby("target_u")["mirna"].nunique()
    g_weight = long_df.groupby("target_u")["mirna_w"].sum()
    idx = g_count.index
    return pd.DataFrame(
        {
            "gene_u": idx.values,
            "n_mirnas": g_count.values,
            "weighted_burden": g_weight.reindex(idx).fillna(0.0).values,
        }
    )


def run_celltype_strata_crossmodal(
    out_dir: Path,
    log,
    *,
    min_genes: int = 120,
    n_perm: int = 500,
    seed: int = 44,
) -> None:
    log("\n=== Cell-type / stratum–resolved cross-modal burden vs aging DE ===")
    man_path = out_dir / "gse188646_strata" / "manifest.csv"
    long_path = out_dir / "mirna_targets_long.csv"
    mirna_sum = out_dir / "mirna_htnsc_astrocyte_summary.csv"
    if not man_path.is_file():
        log(
            "Cell-type cross-modal: no outputs/gse188646_strata/manifest.csv — "
            "run stratified pseudobulk (GSE188646_STRATUM_COL=seurat_clusters or "
            "GSE188646_AUTO_STRATIFIED_PSEUDOBULK=1 with RDS present). Skipped."
        )
        return
    if not long_path.is_file() or not mirna_sum.is_file():
        log("Cell-type cross-modal: missing mirna_targets_long or mirna summary; skipped.")
        return

    mirna_w = pd.read_csv(mirna_sum)
    mirna_w["mirna"] = mirna_w["mirna"].astype(str).str.strip()
    wmap = mirna_w.set_index("mirna")["logfc_htnsc_vs_astro"].to_dict()

    long_df = pd.read_csv(long_path)
    long_df["mirna"] = long_df["mirna"].astype(str).str.strip()
    long_df["target_u"] = long_df["target_gene"].astype(str).str.strip().str.upper()
    genes_burden = _build_gene_burden(long_df, wmap)
    union_set = set(long_df["target_u"])

    manifest = pd.read_csv(man_path)
    if manifest.empty:
        log("Cell-type cross-modal: empty manifest; skipped.")
        return

    rng = np.random.default_rng(seed)
    rows: list[dict] = []

    for _, mr in manifest.iterrows():
        deg_fn = str(mr.get("deg_csv", "")).strip()
        if not deg_fn:
            continue
        deg_path = out_dir / "gse188646_strata" / deg_fn
        if not deg_path.is_file():
            continue
        deg = pd.read_csv(deg_path)
        if "gene" not in deg.columns or "logFC" not in deg.columns:
            continue
        deg["gene_u"] = deg["gene"].astype(str).str.strip().str.upper()
        m = deg.merge(genes_burden, on="gene_u", how="inner")
        if len(m) < min_genes:
            continue

        rho, pval = _spearman(m["weighted_burden"].values, m["logFC"].values)
        y = m["logFC"].values.astype(float)
        w = m["weighted_burden"].values.astype(float)
        null = []
        for _ in range(n_perm):
            yp = rng.permutation(y)
            r0, _ = _spearman(w, yp)
            null.append(r0)
        null = np.asarray(null, dtype=float)
        perm_p = float(np.mean(np.abs(null) >= abs(rho))) if np.isfinite(rho) else float("nan")

        in_union = deg["gene_u"].isin(union_set)
        u_fc = np.abs(deg.loc[in_union, "logFC"].astype(float).values)
        o_fc = np.abs(deg.loc[~in_union, "logFC"].astype(float).values)
        try:
            mw_p = float(stats.mannwhitneyu(u_fc, o_fc, alternative="two-sided").pvalue)
        except ValueError:
            mw_p = float("nan")

        rows.append(
            {
                "stratum_col": str(mr.get("stratum_col", "")),
                "stratum": str(mr.get("stratum", "")),
                "n_cells": int(mr.get("n_cells", 0)) if pd.notna(mr.get("n_cells")) else None,
                "n_young_rep": int(mr.get("n_young_rep", 0)) if pd.notna(mr.get("n_young_rep")) else None,
                "n_aged_rep": int(mr.get("n_aged_rep", 0)) if pd.notna(mr.get("n_aged_rep")) else None,
                "n_genes_merged": int(len(m)),
                "spearman_rho_weighted_burden_vs_logFC": rho,
                "spearman_p_weighted": pval,
                "perm_p_rho_gene_shuffle": perm_p,
                "perm_n": n_perm,
                "mannwhitney_abs_logFC_union_vs_nonunion_p": mw_p,
                "deg_csv": deg_fn,
            }
        )

    if not rows:
        log(
            f"Cell-type cross-modal: no strata passed merge (min_genes={min_genes}). "
            "Relax GSE188646_STRATUM_MIN_CELLS / min_rep in R or lower min_genes in Python."
        )
        return

    out_csv = out_dir / "exploratory_crossmodal_celltype_strata_summary.csv"
    pd.DataFrame(rows).sort_values("spearman_rho_weighted_burden_vs_logFC", ascending=False).to_csv(
        out_csv, index=False
    )
    summ = {
        "n_strata_tested": len(rows),
        "min_genes_per_stratum": min_genes,
        "n_perm_gene_shuffle": n_perm,
        "methodology_note": (
            "Per-stratum pseudobulk DE (young vs aged) within Seurat metadata stratum; "
            "burden is identical to global cross-modal construction. Gene-shuffle null permutes "
            "logFC labels within the merged gene list for that stratum only."
        ),
        "caveat": (
            "Strata are author clusters or user-supplied metadata columns, not necessarily pure "
            "cell types; sex/assay structure matches parent object."
        ),
    }
    (out_dir / "exploratory_crossmodal_celltype_strata_summary.json").write_text(
        json.dumps(summ, indent=2), encoding="utf-8"
    )
    log(f"Wrote {out_csv.name} ({len(rows)} strata).")
    log(json.dumps(summ, indent=2))
