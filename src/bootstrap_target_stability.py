"""Bootstrap stability of miRNA→target union under resampling of top-ranked miRNAs."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mirna_target_union import union_targets_for_mirnas


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def bootstrap_target_stability(
    mirna_ranked_df: pd.DataFrame,
    mirna_col: str,
    score_col: str,
    gmt_path,
    pool_size: int = 80,
    n_draw: int = 40,
    sample_size: int = 35,
    seed: int = 7,
) -> pd.DataFrame:
    """
    Pool = top `pool_size` miRNAs by score. Each bootstrap: draw `sample_size` miRNAs
    with replacement from pool, build miRTarBase union, Jaccard vs reference union
    from top `sample_size` distinct miRNAs (no replacement) for fair comparison — actually
    use fixed top `sample_size` miRNAs as reference union.
    """
    rng = np.random.default_rng(seed)
    df = mirna_ranked_df.sort_values(score_col, ascending=False).reset_index(drop=True)
    n_avail = len(df)
    if n_avail == 0:
        raise ValueError("mirna_ranked_df is empty after sort")
    pool_size = min(pool_size, n_avail)
    sample_size = min(sample_size, pool_size)
    pool = df.head(pool_size)[mirna_col].astype(str).tolist()
    ref_mir = df.head(sample_size)[mirna_col].astype(str).tolist()
    ref_union, _ = union_targets_for_mirnas(ref_mir, gmt_path)
    ref_set = {g.upper() for g in ref_union}

    jac = []
    sizes = []
    for _ in range(n_draw):
        draw = list(rng.choice(pool, size=sample_size, replace=True))
        u, _ = union_targets_for_mirnas(draw, gmt_path)
        s = {g.upper() for g in u}
        jac.append(jaccard(ref_set, s))
        sizes.append(len(s))
    out = pd.DataFrame(
        {
            "jaccard_vs_top_ref": jac,
            "union_size": sizes,
        }
    )
    out.attrs["summary"] = {
        "mean_jaccard": float(np.mean(jac)),
        "median_jaccard": float(np.median(jac)),
        "q05_jaccard": float(np.quantile(jac, 0.05)),
        "q95_jaccard": float(np.quantile(jac, 0.95)),
        "ref_union_size": len(ref_set),
        "pool_size": pool_size,
        "n_draw": n_draw,
        "sample_size": sample_size,
    }
    return out


def save_bootstrap(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)
    summ = df.attrs.get("summary", {})
    if summ:
        lines = [f"{k}: {v}" for k, v in summ.items()]
        path.with_name(path.stem + "_summary.txt").write_text("\n".join(lines), encoding="utf-8")
