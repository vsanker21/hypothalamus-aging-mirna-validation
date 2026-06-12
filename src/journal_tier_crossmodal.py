"""
Cross-modal gene-level analyses linking htNSC miRNA target structure to GSE188646 aging DE.

Exploratory / hypothesis-generating only: miRNA layer (in vitro NSC) and hypothalamic aging DE
are not exchangeable experiments. Outputs use the exploratory_* prefix where appropriate.

Computes:
  - Per-gene target burden (count and htNSC-logFC-weighted sum of incoming miRNA regulators)
  - Spearman correlation burden vs pseudobulk logFC on gene intersection
  - Gene-label permutation null for Spearman rho
  - Random miRNA-set null: repeat burden construction with random miRNAs matched in count from MOESM table
  - Mann–Whitney test of |logFC| in target-union genes vs non-target genes (same universe)
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from gmt_io import iter_gmt_lines
from mirna_target_union import ensure_mirtarbase_gmt


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    r, p = stats.spearmanr(x, y, nan_policy="omit")
    if not np.isfinite(r):
        return float("nan"), float("nan")
    return float(r), float(p)


def run_crossmodal_mirna_aging(
    out_dir: Path,
    log,
    *,
    n_perm: int = 1000,
    n_mirna_draws: int = 200,
    seed: int = 42,
) -> None:
    log("\n=== Journal-tier cross-modal: miRNA target burden vs GSE188646 aging logFC ===")
    deg_path = out_dir / "gse188646_young_vs_aged_deg.csv"
    long_path = out_dir / "mirna_targets_long.csv"
    mirna_sum = out_dir / "mirna_htnsc_astrocyte_summary.csv"
    if not deg_path.is_file() or not long_path.is_file() or not mirna_sum.is_file():
        log("Cross-modal: missing DE, mirna_targets_long, or mirna summary; skipped.")
        return

    deg = pd.read_csv(deg_path)
    if "gene" not in deg.columns or "logFC" not in deg.columns:
        log("Cross-modal: unexpected DE columns; skipped.")
        return
    deg["gene_u"] = deg["gene"].astype(str).str.strip().str.upper()
    se_col = "se_logFC" if "se_logFC" in deg.columns else None

    mirna_w = pd.read_csv(mirna_sum)
    mirna_w["mirna"] = mirna_w["mirna"].astype(str).str.strip()
    wmap = mirna_w.set_index("mirna")["logfc_htnsc_vs_astro"].to_dict()

    long_df = pd.read_csv(long_path)
    long_df["mirna"] = long_df["mirna"].astype(str).str.strip()
    long_df["target_u"] = long_df["target_gene"].astype(str).str.strip().str.upper()

    # Weight per edge: max(0, htNSC vs astro logFC) for that miRNA
    long_df["mirna_w"] = long_df["mirna"].map(wmap).fillna(0.0).clip(lower=0.0)

    g_count = long_df.groupby("target_u")["mirna"].nunique()
    g_weight = long_df.groupby("target_u")["mirna_w"].sum()
    idx = g_count.index
    genes = pd.DataFrame(
        {
            "gene_u": idx.values,
            "n_mirnas": g_count.values,
            "weighted_burden": g_weight.reindex(idx).fillna(0.0).values,
        }
    )
    m = deg.merge(genes, on="gene_u", how="inner")
    if len(m) < 200:
        log(f"Cross-modal: too few genes after merge ({len(m)}); skipped.")
        return

    rho, pval = _spearman(m["weighted_burden"].values, m["logFC"].values)
    rho_c, p_c = _spearman(m["n_mirnas"].values, m["logFC"].values)

    rng = np.random.default_rng(seed)
    y = m["logFC"].values.astype(float)
    w = m["weighted_burden"].values.astype(float)
    null = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        r0, _ = _spearman(w, yp)
        null.append(r0)
    null = np.array(null, dtype=float)
    perm_p = float(np.mean(np.abs(null) >= abs(rho))) if np.isfinite(rho) else float("nan")

    # |logFC| in union targets vs other genes in DE table (not necessarily in GMT universe)
    in_union = deg["gene_u"].isin(set(long_df["target_u"]))
    u_fc = np.abs(deg.loc[in_union, "logFC"].astype(float).values)
    o_fc = np.abs(deg.loc[~in_union, "logFC"].astype(float).values)
    try:
        mw = stats.mannwhitneyu(u_fc, o_fc, alternative="two-sided")
        mw_p = float(mw.pvalue)
        mw_p_greater = float(stats.mannwhitneyu(u_fc, o_fc, alternative="greater").pvalue)
        mw_p_less = float(stats.mannwhitneyu(u_fc, o_fc, alternative="less").pvalue)
    except ValueError:
        mw_p = mw_p_greater = mw_p_less = float("nan")
    med_u, med_o = float(np.median(u_fc)), float(np.median(o_fc))
    mean_u, mean_o = float(np.mean(u_fc)), float(np.mean(o_fc))

    # Random miRNA sets (same cardinality as distinct miRNAs in long_df)
    mirnas_obs = sorted(long_df["mirna"].unique())
    n_mir = len(mirnas_obs)
    gmt = ensure_mirtarbase_gmt()
    full_mir_targets: dict[str, set[str]] = {}
    for term, genes in iter_gmt_lines(gmt):
        if not (term.startswith("mmu-miR") or term.startswith("mmu-let")):
            continue
        if term not in wmap:
            continue
        full_mir_targets[term] = {str(g).strip().upper() for g in genes if str(g).strip()}
    pool_draw = sorted(full_mir_targets.keys())
    if len(pool_draw) < n_mir + 5:
        log("Cross-modal: insufficient GMT×MOESM miRNA pool for random-set null; skipping draws.")
        rand_rhos = np.array([], dtype=float)
    else:
        rand_rhos = []
        for _ in range(n_mirna_draws):
            draw = list(rng.choice(pool_draw, size=n_mir, replace=False))
            gene_w: dict[str, float] = defaultdict(float)
            for mir in draw:
                ww = float(max(0.0, wmap.get(mir, 0.0)))
                for g in full_mir_targets.get(mir, ()):
                    gene_w[g] += ww
            if not gene_w:
                continue
            tmp = pd.DataFrame([{"gene_u": g, "wb_rand": v} for g, v in gene_w.items()])
            mm = deg.merge(tmp, on="gene_u", how="inner")
            if len(mm) < 200:
                continue
            rr, _ = _spearman(mm["wb_rand"].values, mm["logFC"].values)
            if np.isfinite(rr):
                rand_rhos.append(rr)
        rand_rhos = np.asarray(rand_rhos, dtype=float)
    mirna_null_p = float(np.mean(np.abs(rand_rhos) >= abs(rho))) if rand_rhos.size and np.isfinite(rho) else float("nan")

    out_gene = m[["gene_u", "logFC", "n_mirnas", "weighted_burden"]].rename(columns={"gene_u": "gene"})
    if se_col:
        out_gene["se_logFC"] = m[se_col].values
    out_gene.to_csv(out_dir / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv", index=False)

    summ = {
        "n_genes_merged": int(len(m)),
        "n_distinct_mirnas_in_program": int(n_mir),
        "spearman_rho_weighted_burden_vs_logFC": rho,
        "spearman_p_weighted": pval,
        "spearman_rho_count_vs_logFC": rho_c,
        "spearman_p_count": p_c,
        "perm_p_rho_gene_shuffle_n": n_perm,
        "perm_p_rho_gene_shuffle": perm_p,
        "mannwhitney_abs_logFC_union_vs_nonunion_p": mw_p,
        "mannwhitney_abs_logFC_union_greater_p": mw_p_greater,
        "mannwhitney_abs_logFC_union_less_p": mw_p_less,
        "median_abs_logFC_union": med_u,
        "median_abs_logFC_nonunion": med_o,
        "mean_abs_logFC_union": mean_u,
        "mean_abs_logFC_nonunion": mean_o,
        "delta_median_abs_logFC_union_minus_nonunion": med_u - med_o,
        "random_mirna_set_draws": int(n_mirna_draws),
        "random_mirna_set_abs_rho_ge_observed_frac": mirna_null_p,
        "caveat": (
            "MOESM htNSC program vs hypothalamus aging DE; not causal; exploratory only. "
            "Two-sided Mann–Whitney on |logFC| tests distributional difference; in GSE188646 pseudobulk "
            "union targets have lower median/mean |logFC| (attenuated volatility), not larger shifts."
        ),
    }
    (out_dir / "exploratory_crossmodal_mirna_aging_summary.json").write_text(
        json.dumps(summ, indent=2), encoding="utf-8"
    )
    pd.DataFrame({"perm_rho_weighted_burden": null}).to_csv(
        out_dir / "exploratory_crossmodal_permutation_rho_gene_shuffle.csv", index=False
    )
    if rand_rhos.size:
        pd.DataFrame({"random_set_rho": rand_rhos}).to_csv(
            out_dir / "exploratory_crossmodal_random_mirna_set_rho_null.csv", index=False
        )

    log(json.dumps(summ, indent=2))
