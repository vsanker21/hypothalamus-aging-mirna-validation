"""
Second-pass negative controls for burden vs aging logFC (Science Advances-style rigor).

Complements journal_tier_crossmodal.py:
  - Unrestricted gene-label shuffle ignores confounding by detection / hub targeting.
  - Here we shuffle logFC only within strata matched on:
      (1) DE precision (se_logFC deciles),
      (2) miRNA-program target degree (n_mirnas, capped),
      (3) combined degree × precision,
      (4) GMT-wide weighted in-degree from the same MOESM-weighted miRNA pool (all mmu terms in GMT
          with a MOESM weight) × precision — a "network exposure" control for how targetable a gene is
          across the whole miRTarBase layer, not only the htNSC top program.

RATIONALE (bipartite nulls vs this implementation)
  A *strict* degree-preserving null on the full miRNA–gene bipartite graph would sample (or MCMC
  swap) edges while fixing each miRNA's out-degree sequence and each gene's in-degree sequence
  (bipartite configuration model / stub matching). That is the gold standard when edges are
  treated as exchangeable *and* the graph is small enough for many exact draws. Here, edges are
  literature-supported priors (not exchangeable draws), the graph is large, and we keep the
  *observed* program burden vector w fixed to test its alignment with aging logFC.

  We therefore use the usual omics-scale compromise: **coarsen** continuous summaries (DE
  precision via se_logFC deciles; program in-degree capped; GMT-wide weighted in-degree deciles)
  and apply **restricted permutations** of logFC within strata. This destroys gene-level coupling
  between y and w while preserving marginal structure that would otherwise inflate false positives
  when highly targetable or precisely estimated genes also drift in logFC.

  GMT-wide indegree adds a **network-matched exposure** dimension beyond the htNSC top program,
  separate from random miRNA-*set* resampling in journal_tier_crossmodal.py.

Exploratory only; same cross-study caveats as crossmodal.
"""
from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from gmt_io import iter_gmt_lines
from mirna_target_union import ensure_mirtarbase_gmt


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    r, _ = stats.spearmanr(x, y, nan_policy="omit")
    return float(r) if np.isfinite(r) else float("nan")


def _stratified_shuffle_y(
    y: np.ndarray,
    strata: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Within each stratum id, shuffle y; leave singleton strata unchanged."""
    out = y.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        if len(idx) <= 1:
            continue
        out[idx] = rng.permutation(out[idx])
    return out


def _make_se_bins(se: np.ndarray, q: int = 10) -> np.ndarray:
    s = pd.Series(se)
    try:
        return pd.qcut(s.rank(method="first"), q=q, labels=False, duplicates="drop").astype(int).values
    except ValueError:
        return pd.cut(s, bins=min(q, max(3, int(s.nunique()))), labels=False, duplicates="drop").astype(int).values


def _gmt_weighted_indegree(gmt_path: Path, wmap: dict[str, float]) -> dict[str, float]:
    """Sum of MOESM weights for every mmu miRNA in GMT that lists the gene."""
    indeg: dict[str, float] = defaultdict(float)
    for term, genes in iter_gmt_lines(gmt_path):
        if not (term.startswith("mmu-miR") or term.startswith("mmu-let")):
            continue
        ww = float(max(0.0, wmap.get(term, 0.0)))
        if ww <= 0:
            continue
        for g in genes:
            indeg[str(g).strip().upper()] += ww
    return indeg


def run_expression_degree_negative_controls(
    out_dir: Path,
    log,
    *,
    n_perm: int = 800,
    seed: int = 43,
) -> None:
    log("\n=== Journal-tier negative controls: precision (se_logFC) + degree-matched logFC shuffles ===")
    deg_path = out_dir / "gse188646_young_vs_aged_deg.csv"
    long_path = out_dir / "mirna_targets_long.csv"
    mirna_sum = out_dir / "mirna_htnsc_astrocyte_summary.csv"
    if not deg_path.is_file() or not long_path.is_file() or not mirna_sum.is_file():
        log("Negative controls: missing DE, mirna_targets_long, or mirna summary; skipped.")
        return

    deg = pd.read_csv(deg_path)
    if "gene" not in deg.columns or "logFC" not in deg.columns or "se_logFC" not in deg.columns:
        log("Negative controls: need se_logFC in DE table for precision strata; skipped.")
        return

    deg["gene_u"] = deg["gene"].astype(str).str.strip().str.upper()
    mirna_w = pd.read_csv(mirna_sum)
    mirna_w["mirna"] = mirna_w["mirna"].astype(str).str.strip()
    wmap = mirna_w.set_index("mirna")["logfc_htnsc_vs_astro"].to_dict()

    long_df = pd.read_csv(long_path)
    long_df["mirna"] = long_df["mirna"].astype(str).str.strip()
    long_df["target_u"] = long_df["target_gene"].astype(str).str.strip().str.upper()
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
        log(f"Negative controls: too few genes after merge ({len(m)}); skipped.")
        return

    w = m["weighted_burden"].values.astype(float)
    y_obs = m["logFC"].values.astype(float)
    se = m["se_logFC"].values.astype(float)
    n_mir = m["n_mirnas"].values.astype(int)
    rho_obs = _spearman(w, y_obs)

    gmt = ensure_mirtarbase_gmt()
    gmt_w = _gmt_weighted_indegree(gmt, wmap)
    m["gmt_weighted_indegree"] = m["gene_u"].map(lambda g: float(gmt_w.get(g, 0.0)))

    rng = np.random.default_rng(seed)

    se_bin = _make_se_bins(se, q=10)
    deg_cap = np.minimum(n_mir, 12)
    gmt_deg_bin = _make_se_bins(m["gmt_weighted_indegree"].values, q=10)

    # Strata as integer keys (combine with large primes to reduce collision)
    strata_se = se_bin
    strata_deg_se = deg_cap * 1000 + se_bin
    strata_gmt_se = gmt_deg_bin * 1000 + se_bin

    rows = []
    null_se = []
    null_deg_se = []
    null_gmt_se = []
    for it in range(n_perm):
        y1 = _stratified_shuffle_y(y_obs, strata_se, rng)
        y2 = _stratified_shuffle_y(y_obs, strata_deg_se, rng)
        y3 = _stratified_shuffle_y(y_obs, strata_gmt_se, rng)
        r1, r2, r3 = _spearman(w, y1), _spearman(w, y2), _spearman(w, y3)
        null_se.append(r1)
        null_deg_se.append(r2)
        null_gmt_se.append(r3)
        rows.append({"iter": it, "rho_perm_se_strata": r1, "rho_perm_deg_se_strata": r2, "rho_perm_gmtdeg_se_strata": r3})

    null_se = np.asarray(null_se, dtype=float)
    null_deg_se = np.asarray(null_deg_se, dtype=float)
    null_gmt_se = np.asarray(null_gmt_se, dtype=float)

    def p_two(rnull: np.ndarray, robs: float) -> float:
        if not np.isfinite(robs):
            return float("nan")
        rnull = rnull[np.isfinite(rnull)]
        if rnull.size == 0:
            return float("nan")
        return float(np.mean(np.abs(rnull) >= abs(robs)))

    summ = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "n_genes": int(len(m)),
        "spearman_rho_observed_weighted_burden_vs_logFC": rho_obs,
        "n_perm": int(n_perm),
        "perm_p_abs_rho_ge_obs_se_decile_strata_only": p_two(null_se, rho_obs),
        "perm_p_abs_rho_ge_obs_program_degree_cap_x_se_strata": p_two(null_deg_se, rho_obs),
        "perm_p_abs_rho_ge_obs_gmt_weighted_indegree_x_se_strata": p_two(null_gmt_se, rho_obs),
        "strata_notes": (
            "se_strata: qcut ranks of se_logFC into deciles. "
            "deg_se: capped program n_mirnas (12+) crossed with se decile. "
            "gmt_se: GMT-wide MOESM-weighted targetability decile x se decile."
        ),
        "methodology_note": (
            "Restricted permutation of logFC within coarse strata (precision, program degree, "
            "GMT-wide targetability) approximates adjustment for hub targeting and DE precision "
            "without full bipartite configuration sampling; see module docstring and "
            "data/provenance/JOURNAL_TIER_COMPUTATIONAL_README.txt."
        ),
        "caveat": "Stratified shuffles preserve marginal coupling of logFC to precision/degree bins; "
        "they do not remove all hidden confounders. Exploratory only.",
    }

    pd.DataFrame(rows).to_csv(
        out_dir / "exploratory_negative_controls_stratified_perm_rhos.csv", index=False
    )
    (out_dir / "exploratory_negative_controls_summary.json").write_text(
        json.dumps(summ, indent=2), encoding="utf-8"
    )

    m_out = m[
        ["gene_u", "logFC", "se_logFC", "n_mirnas", "weighted_burden", "gmt_weighted_indegree"]
    ].rename(columns={"gene_u": "gene"})
    m_out.to_csv(out_dir / "exploratory_negative_controls_gene_covariates.csv", index=False)

    log(json.dumps(summ, indent=2))
