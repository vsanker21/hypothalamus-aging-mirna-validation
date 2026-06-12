"""Fisher overlap tests between miRNA target union and external gene lists."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from mirna_target_union import mmu_gmt_gene_universe


def _read_gene_column(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    lower = {c.lower(): c for c in df.columns}
    for key in ("gene", "gene_symbol", "symbol", "external_gene_name"):
        if key in lower:
            col = lower[key]
            out = df[[col]].copy()
            out.columns = ["gene"]
            return out.dropna()
    raise ValueError(f"No gene column found in {path}; expected gene/symbol/external_gene_name")


def _read_seurat_markers(path: Path, padj_max: float = 0.05) -> set[str]:
    """Seurat FindMarkers output: row names = genes OR column gene."""
    df = pd.read_csv(path)
    if "gene" not in df.columns and df.index.name != "gene":
        if df.select_dtypes(include=["object"]).shape[1] == 0:
            df = df.reset_index().rename(columns={"index": "gene"})
    if "gene" not in df.columns:
        first = df.columns[0]
        df = df.rename(columns={first: "gene"})
    padj_col = None
    for c in ("p_val_adj", "padj", "FDR", "adj.P.Val"):
        if c in df.columns:
            padj_col = c
            break
    if padj_col:
        df = df[df[padj_col] <= padj_max]
    return set(df["gene"].astype(str).str.strip())


def deg_top_by_rank_score(path: Path, top_n: int) -> set[str]:
    """
    Exploratory gene set: top `top_n` genes by sign(logFC) * (-log10(p)) from a DE table
    with columns gene, logFC, p_val (edgeR/Seurat style). No FDR filter.
    """
    df = pd.read_csv(path)
    if "gene" not in df.columns:
        first = df.columns[0]
        df = df.rename(columns={first: "gene"})
    p_col = "p_val" if "p_val" in df.columns else None
    if p_col is None:
        for c in ("PValue", "pvalue", "P.Value"):
            if c in df.columns:
                p_col = c
                break
    if p_col is None:
        raise ValueError(f"No p-value column in {path}")
    if "logFC" not in df.columns:
        raise ValueError(f"No logFC column in {path}")
    pv = df[p_col].astype(float).clip(lower=1e-300)
    score = np.sign(df["logFC"].astype(float)) * (-np.log10(pv))
    df = df.assign(_s=score).sort_values("_s", ascending=False)
    return set(df.head(int(top_n))["gene"].astype(str).str.strip())


def gse188646_de_calibration(path: Path) -> pd.DataFrame:
    """Counts of genes meeting nominal FDR cutoffs (edgeR column p_val_adj)."""
    df = pd.read_csv(path)
    padj_col = None
    for c in ("p_val_adj", "padj", "FDR", "adj.P.Val"):
        if c in df.columns:
            padj_col = c
            break
    if padj_col is None:
        return pd.DataFrame({"note": ["no padj column"]})
    pv = df[padj_col].astype(float)
    rows = []
    for t in (0.05, 0.1, 0.2, 0.5, 1.0):
        rows.append({"padj_max": t, "n_genes": int((pv <= t).sum())})
    return pd.DataFrame(rows)


def fisher_overlap(
    targets: set[str],
    external: set[str],
    universe: set[str],
) -> dict:
    """
    Fisher exact test on 2x2 contingency table restricted to a common universe U
    (typically miRTarBase mmu target gene background). All symbols compared uppercase.
    T' = T∩U, E' = E∩U; a=|T'∩E'|, b=|T'\\E'|, c=|E'\\T'|, d=|U \\ (T'∪E')|.
    """
    def norm(s: set[str]) -> set[str]:
        return {x.strip().upper() for x in s if x and str(x).strip()}

    universe = norm(universe)
    targets = norm(targets) & universe
    external = norm(external) & universe
    a = len(targets & external)
    b = len(targets - external)
    c = len(external - targets)
    d = len(universe - targets - external)
    if min(a, b, c, d) < 0:
        raise ValueError("negative contingency")
    oddsr, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    or_val = float(oddsr) if np.isfinite(oddsr) else None
    if or_val == 0.0 and min(b, c) > 0:
        or_val = (a * d) / (b * c)
    return {
        "a_targets_and_external": a,
        "b_targets_only": b,
        "c_external_only": c,
        "d_neither": d,
        "universe_size": len(universe),
        "odds_ratio": or_val,
        "fisher_two_sided_p": float(p),
    }


def load_hypomap_de_union(
    urls: list[tuple[str, str]], padj_max: float = 0.05, abs_lfc_min: float = 0.5
) -> tuple[set[str], pd.DataFrame]:
    """Download HypoMap GEO supplementary DESeq2 tables (NOT young-vs-aged)."""
    import gzip
    import io

    import requests

    rows = []
    genes: set[str] = set()
    for label, url in urls:
        r = requests.get(url, timeout=180)
        r.raise_for_status()
        buf = io.BytesIO(r.content)
        with gzip.open(buf, "rt", encoding="utf-8", errors="replace") as f:
            sub = pd.read_csv(f)
        sub = sub.assign(source_table=label)
        if "padj" in sub.columns:
            sub = sub[sub["padj"] <= padj_max]
        if "log2FoldChange" in sub.columns:
            sub = sub[sub["log2FoldChange"].abs() >= abs_lfc_min]
        if "external_gene_name" in sub.columns:
            gcol = "external_gene_name"
        elif "gene" in sub.columns:
            gcol = "gene"
        else:
            raise ValueError(f"No gene column in {label}")
        for g in sub[gcol].dropna().astype(str).unique():
            genes.add(g.strip())
        rows.append(sub)
    meta = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return genes, meta


def save_overlap_report(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
