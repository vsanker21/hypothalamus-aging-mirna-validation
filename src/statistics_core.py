"""Statistical summaries for supplementary tables (small-n aware)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def paired_ttest(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) < 2:
        return {"n": len(a), "mean_diff": np.nan, "t": np.nan, "p_two_sided": np.nan}
    d = a - b
    t, p = stats.ttest_1samp(d, popmean=0.0, alternative="two-sided")
    return {
        "n": int(len(a)),
        "mean_diff": float(np.mean(d)),
        "sd_diff": float(np.std(d, ddof=1)),
        "t": float(t),
        "p_two_sided": float(p),
    }


def one_way_anova_three_groups(
    g1: np.ndarray, g2: np.ndarray, g3: np.ndarray
) -> dict:
    """One-way ANOVA; omits nan per group."""
    g1 = np.asarray(g1, dtype=float)[np.isfinite(g1)]
    g2 = np.asarray(g2, dtype=float)[np.isfinite(g2)]
    g3 = np.asarray(g3, dtype=float)[np.isfinite(g3)]
    if min(len(g1), len(g2), len(g3)) < 2:
        return {"f": np.nan, "p": np.nan}
    f, p = stats.f_oneway(g1, g2, g3)
    return {"f": float(f), "p": float(p), "n1": len(g1), "n2": len(g2), "n3": len(g3)}


def cohens_d_paired(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    d = a[mask] - b[mask]
    if len(d) < 2 or np.std(d, ddof=1) == 0:
        return float("nan")
    return float(np.mean(d) / np.std(d, ddof=1))


def mirna_specificity_null(
    htnsc1: np.ndarray,
    htnsc2: np.ndarray,
    a1: np.ndarray,
    a2: np.ndarray,
    n_perm: int = 8000,
    seed: int = 42,
) -> dict:
    """
    Per-miRNA null: randomly assign the four replicate values into two groups of two
    (labelled 'htNSC' vs 'astrocyte'); take median |mean(H)-mean(A)| across miRNAs.
    Observed statistic uses the true htNSC vs astrocyte pairing.
    """
    rng = np.random.default_rng(seed)
    h1 = np.asarray(htnsc1, dtype=float)
    h2 = np.asarray(htnsc2, dtype=float)
    x1 = np.asarray(a1, dtype=float)
    x2 = np.asarray(a2, dtype=float)
    mask = np.isfinite(h1) & np.isfinite(h2) & np.isfinite(x1) & np.isfinite(x2)
    h1, h2, x1, x2 = h1[mask], h2[mask], x1[mask], x2[mask]
    obs_vals = np.stack([h1, h2, x1, x2], axis=1)
    obs_stat = float(
        np.median(np.abs((h1 + h2) / 2.0 - (x1 + x2) / 2.0))
    )
    perm_stats = np.empty(n_perm)
    mat = obs_vals.copy()
    for i in range(n_perm):
        sh = rng.permuted(mat, axis=1)
        g1 = sh[:, :2].mean(axis=1)
        g2 = sh[:, 2:].mean(axis=1)
        perm_stats[i] = float(np.median(np.abs(g1 - g2)))
    p = float(np.mean(perm_stats >= obs_stat))
    return {
        "median_abs_logfc_htnsc_vs_astro": obs_stat,
        "perm_p_median_vs_null": p,
        "n_mirna": int(len(obs_vals)),
        "n_perm": n_perm,
    }
