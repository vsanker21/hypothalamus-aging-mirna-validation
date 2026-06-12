"""Build union of experimentally supported mRNA targets for top htNSC-enriched miRNAs."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from config import DATA_DIR, OUTPUT_DIR
from gmt_io import CACHE_DIR, build_mirna_to_targets_mmu, download_text, iter_gmt_lines

MIRTARBASE_GMT_URL = (
    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=miRTarBase_2017"
)
GMT_CACHE = CACHE_DIR / "miRTarBase_2017_Enrichr.gmt"


def load_top_htnsc_mirnas(
    summary_csv: Path, top_n: int = 60, min_logfc: float = 0.5
) -> list[str]:
    df = pd.read_csv(summary_csv)
    df = df.sort_values("logfc_htnsc_vs_astro", ascending=False)
    df = df[df["logfc_htnsc_vs_astro"] >= min_logfc]
    top = df.head(top_n)["mirna"].astype(str).tolist()
    return top


def ensure_mirtarbase_gmt() -> Path:
    if not GMT_CACHE.exists():
        download_text(MIRTARBASE_GMT_URL, GMT_CACHE, timeout=600)
    return GMT_CACHE


def union_targets_for_mirnas(
    mirnas: list[str], gmt_path: Path | None = None
) -> tuple[set[str], dict[str, set[str]]]:
    gmt_path = gmt_path or ensure_mirtarbase_gmt()
    allow = set(mirnas)
    mt = build_mirna_to_targets_mmu(gmt_path, lambda t: t in allow)
    union: set[str] = set()
    for gset in mt.values():
        union |= gset
    return union, mt


def mmu_gmt_gene_universe(gmt_path: Path) -> set[str]:
    """All genes appearing under any mmu-miR / mmu-let term in miRTarBase GMT (Enrichr)."""
    u: set[str] = set()
    for term, genes in iter_gmt_lines(gmt_path):
        if term.startswith("mmu-miR") or term.startswith("mmu-let"):
            u.update(genes)
    return u


def save_target_maps(
    union: set[str], per_mirna: dict[str, set[str]], out_dir: Path
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.Series(sorted(union), name="gene").to_csv(
        out_dir / "mirna_target_union_genes.csv", index=False
    )
    rows = []
    for mir, genes in sorted(per_mirna.items()):
        for g in sorted(genes):
            rows.append({"mirna": mir, "target_gene": g})
    pd.DataFrame(rows).to_csv(out_dir / "mirna_targets_long.csv", index=False)


def build_ranked_genes_for_gsea(
    per_mirna: dict[str, set[str]], universe: set[str]
) -> pd.DataFrame:
    """
    Rank genes by number of top-miRNA regulators (restricted to miRTarBase mmu universe).
    Genes never targeted get rank_metric 0 so the preranked list spans the full universe.
    Tie-break with tiny noise so gseapy prerank does not see massive duplicate ranks.
    """
    import numpy as np

    rng = np.random.default_rng(42)
    counts: dict[str, float] = {g: 0.0 for g in universe}
    for genes in per_mirna.values():
        for g in genes:
            if g in counts:
                counts[g] += 1.0
    df = pd.DataFrame([{"gene": g, "rank_metric": c} for g, c in counts.items()])
    df = df.sort_values("rank_metric", ascending=False).reset_index(drop=True)
    df["rank_metric"] = df["rank_metric"].values + 1e-9 * rng.random(len(df))
    return df
