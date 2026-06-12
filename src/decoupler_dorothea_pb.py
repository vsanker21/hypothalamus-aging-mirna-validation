"""DoRothEA TF activities on GSE188646 pseudobulk via decoupler (mlm; exploratory)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


def _load_pseudobulk_matrix(
    out_dir: Path,
) -> tuple[sparse.csr_matrix, list[str], list[str]] | None:
    mtx = out_dir / "gse188646_pseudobulk_counts.mtx"
    rows = out_dir / "gse188646_pseudobulk_counts_rownames.csv"
    cols = out_dir / "gse188646_pseudobulk_counts_colnames.csv"
    if not mtx.is_file() or not rows.is_file() or not cols.is_file():
        return None
    X = mmread(mtx)
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()
    g = pd.read_csv(rows)["gene"].astype(str).tolist()
    s = pd.read_csv(cols)["sample"].astype(str).tolist()
    if X.shape[0] != len(g) or X.shape[1] != len(s):
        raise ValueError(
            f"Pseudobulk matrix shape {X.shape} mismatches gene ({len(g)}) or sample ({len(s)}) lists"
        )
    return X, g, s


def run_dorothea_on_pseudobulk(
    out_dir: Path,
    metadata_csv: Path | None = None,
    confidence_levels: tuple[str, ...] = ("A", "B"),
    tmin: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame | None] | None:
    """
    Same log-CPM sample x gene matrix as PROGENy; decoupler mlm on DoRothEA edges.
    (mlm-only: consensus across mlm+ulm can fail scipy FDR when padj contain NaNs on small N.)
    Returns (wide table: sample + per-TF scores, young_vs_aged tests) or None.
    """
    import decoupler as dc

    loaded = _load_pseudobulk_matrix(out_dir)
    if loaded is None:
        return None
    X, genes, samples = loaded
    meta_path = metadata_csv or (out_dir / "gse188646_pseudobulk_metadata.csv")
    if not meta_path.is_file():
        return None
    meta = pd.read_csv(meta_path)
    if "orig.ident" not in meta.columns or "age_bin" not in meta.columns:
        return None
    meta = meta.drop_duplicates(subset=["orig.ident"]).set_index("orig.ident")
    order = [str(x) for x in samples]
    meta = meta.reindex(order)
    if meta["age_bin"].isna().any():
        return None

    xm = np.asarray(X.toarray(), dtype=float)
    libs = xm.sum(axis=0)
    libs = np.where(libs <= 0, np.nan, libs)
    cpm = (xm / libs) * 1e6
    cpm = np.nan_to_num(cpm, nan=0.0, posinf=0.0, neginf=0.0)
    logcpm = np.log2(cpm + 1.0)
    mat = pd.DataFrame(logcpm.T, index=order, columns=[str(g) for g in genes])

    net = dc.op.dorothea(
        organism="mouse",
        levels=list(confidence_levels),
        license="academic",
        verbose=False,
    )
    acts, _pvals = dc.mt.mlm(mat, net, tmin=tmin, verbose=False)
    acts = acts.copy()
    acts.insert(0, "sample", acts.index.astype(str))

    tests = None
    m_y = meta["age_bin"].astype(str).str.strip().eq("Young")
    m_a = meta["age_bin"].astype(str).str.strip().eq("Aged")
    if m_y.sum() >= 2 and m_a.sum() >= 2:
        from scipy.stats import ttest_ind
        from statsmodels.stats.multitest import multipletests

        score_df = acts.set_index("sample")
        m_y = meta["age_bin"].astype(str).str.strip().eq("Young")
        m_a = meta["age_bin"].astype(str).str.strip().eq("Aged")
        rows = []
        for tf in score_df.columns:
            v_y = score_df.loc[m_y, tf].astype(float)
            v_a = score_df.loc[m_a, tf].astype(float)
            tstat, p = ttest_ind(v_y, v_a, equal_var=False)
            rows.append(
                {
                    "tf": tf,
                    "mean_young": float(v_y.mean()),
                    "mean_aged": float(v_a.mean()),
                    "welch_t": float(tstat),
                    "p_two_sided": float(p),
                }
            )
        tests = pd.DataFrame(rows)
        if len(tests):
            _, q, _, _ = multipletests(tests["p_two_sided"].astype(float), method="fdr_bh")
            tests["fdr_bh"] = q

    return acts, tests
