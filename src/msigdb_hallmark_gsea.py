"""MSigDB Hallmark ORA + preranked GSEA (local GMT from Enrichr); aging-signature ORA."""
from __future__ import annotations

import time
from pathlib import Path

import gseapy as gp
import pandas as pd

from gmt_io import CACHE_DIR, download_text

HALLMARK_GMT_URL = (
    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
)
HALLMARK_GMT_CACHE = CACHE_DIR / "MSigDB_Hallmark_2020_Enrichr.gmt"


def ensure_hallmark_gmt() -> Path:
    if not HALLMARK_GMT_CACHE.exists():
        download_text(HALLMARK_GMT_URL, HALLMARK_GMT_CACHE, timeout=300)
    return HALLMARK_GMT_CACHE


def _align_local_ora_like_enrichr(df: pd.DataFrame) -> pd.DataFrame:
    """Add optional Enrichr columns so downstream CSVs stay schema-stable."""
    out = df.copy()
    for c in ("Old P-value", "Old Adjusted P-value"):
        if c not in out.columns:
            out[c] = float("nan")
    return out


def run_hallmark_ora(
    genes: list[str],
    label: str,
    *,
    background_genes: set[str] | None = None,
    max_retries: int = 5,
) -> pd.DataFrame:
    """
    Hallmark ORA via Enrichr API when available; otherwise local hypergeometric ORA
    on the cached MSigDB Hallmark GMT (same file as prerank GSEA).

    When Enrichr is down, ``background_genes`` (typically miRTarBase mmu universe)
    is passed to gseapy.enrich so the null matches overlap/Fisher-style gene space;
    if omitted, gseapy uses the union of genes in the GMT (less ideal).
    """
    genes_clean = sorted({str(g).strip().upper() for g in genes if str(g).strip()})
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            enr = gp.enrichr(
                gene_list=genes_clean,
                gene_sets=["MSigDB_Hallmark_2020"],
                organism="mouse",
                outdir=None,
                no_plot=True,
            )
            df = enr.results.copy()
            if _enrichr_results_look_valid(df):
                df.insert(0, "analysis", label)
                df["ora_backend"] = "enrichr_api"
                return df
            raise RuntimeError("Enrichr hallmark ORA returned invalid table")
        except Exception as e:
            last_err = e
            if attempt + 1 < max_retries:
                time.sleep(min(30.0, 2.0 ** attempt))
    gmt_path = ensure_hallmark_gmt()
    bg: list[str] | None = None
    if background_genes is not None and len(background_genes) > 200:
        bg = sorted({str(g).strip().upper() for g in background_genes})
    try:
        enr2 = gp.enrich(
            gene_list=genes_clean,
            gene_sets=str(gmt_path),
            background=bg,
            no_plot=True,
            verbose=False,
        )
    except Exception as e:
        raise RuntimeError(
            f"Hallmark ORA: Enrichr failed ({last_err!r}) and local gseapy.enrich failed ({e!r})"
        ) from e
    df2 = enr2.results.copy()
    if df2 is None or len(df2) == 0:
        raise RuntimeError("Local Hallmark ORA produced empty results") from last_err
    df2 = _align_local_ora_like_enrichr(df2)
    df2.insert(0, "analysis", label)
    df2["ora_backend"] = "local_gmt_hypergeom"
    df2["ora_background_note"] = (
        f"miRTarBase_mmu_universe_n={len(bg)}" if bg is not None else "gmt_gene_union_default"
    )
    return df2


def _enrichr_results_look_valid(df: pd.DataFrame) -> bool:
    """Reject HTML error pages / empty frames from Enrichr HTTP failures."""
    if df is None or len(df) == 0:
        return False
    col = "Adjusted P-value" if "Adjusted P-value" in df.columns else None
    if col is None:
        return False
    s = df[col].astype(str)
    if s.str.contains("doctype|html|NullPointer|Exception", case=False, na=False).any():
        return False
    return bool(pd.to_numeric(df[col], errors="coerce").notna().any())


def _aging_ora_from_multi_library_terms(terms_path: Path, label: str) -> pd.DataFrame | None:
    """
    Offline fallback: same gene set was already run in multi-library ORA; reuse
    Aging_Perturbations_from_GEO_up rows (identical GMT) when standalone Enrichr fails.
    """
    if not terms_path.is_file():
        return None
    t = pd.read_csv(terms_path)
    if "library" not in t.columns:
        return None
    sub = t.loc[t["library"].astype(str) == "Aging_Perturbations_from_GEO_up"].copy()
    if len(sub) == 0:
        return None
    sub.insert(0, "analysis", label)
    return sub


def run_aging_geo_up_ora(
    genes: list[str],
    label: str,
    *,
    multi_library_terms: Path | None = None,
    max_retries: int = 5,
) -> pd.DataFrame:
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=["Aging_Perturbations_from_GEO_up"],
                organism="mouse",
                outdir=None,
                no_plot=True,
            )
            df = enr.results.copy()
            if _enrichr_results_look_valid(df):
                df.insert(0, "analysis", label)
                return df
            raise RuntimeError("Enrichr returned non-tabular or invalid Adjusted P-value")
        except Exception as e:
            last_err = e
            if attempt + 1 < max_retries:
                time.sleep(min(30.0, 2.0 ** attempt))
    if multi_library_terms is not None:
        fb = _aging_ora_from_multi_library_terms(multi_library_terms, label)
        if fb is not None and len(fb):
            return fb
    if last_err is not None:
        raise last_err
    raise RuntimeError("Aging_Perturbations_from_GEO_up ORA failed with no exception detail")


def run_prerank_hallmark_gsea(
    ranked: pd.DataFrame,
    gene_col: str,
    score_col: str,
    out_stem: Path,
) -> pd.DataFrame | None:
    """
    Preranked GSEA vs MSigDB Hallmark (GMT from Enrichr; gene symbols largely MGI-compatible).
    ranked: two columns — gene symbol and a continuous ranking metric (higher = stronger).
    """
    gmt = ensure_hallmark_gmt()
    rnk = ranked[[gene_col, score_col]].dropna().drop_duplicates(subset=[gene_col])
    if len(rnk) < 30:
        return None
    rnk = rnk.set_index(gene_col)
    rnk.columns = ["rank"]
    outdir = out_stem.parent / f"{out_stem.name}_gsea_prerank"
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=str(gmt),
        outdir=str(outdir),
        permutation_num=999,
        seed=42,
        verbose=False,
        no_plot=True,
        min_size=10,
        max_size=800,
    )
    return pre_res.res2d


def save_df(df: pd.DataFrame | None, path: Path) -> None:
    if df is None or len(df) == 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
