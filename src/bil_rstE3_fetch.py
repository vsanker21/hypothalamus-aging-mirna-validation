"""
Download/cache Jin et al. 2025 RSTE3 processed tables from Brain Image Library (doi:10.35077/g.1157).

Cache root: data/references/bil_rstE3/
"""
from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import requests

BIL_DOI = "10.35077/g.1157"
BIL_DOWNLOAD_BASE = "https://download.brainimagelibrary.org/group/20240509/RSTE3/"

RSTE3_FILES: dict[str, str] = {
    "metadata": "RSTE3_metadata.csv",
    "cellxgene": "RSTE3_cellxgene.csv",
}


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def cache_dir() -> Path:
    return _project_root() / "data" / "references" / "bil_rstE3"


def _download_file(url: str, dest: Path, *, log=None, chunk_mb: int = 1) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 0:
        if log:
            log(f"  cached {dest.name} ({dest.stat().st_size // (1 << 20)} MB)")
        return dest

    if log:
        log(f"  downloading {dest.name} from BIL...")
    t0 = time.time()
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_mb << 20):
                if chunk:
                    f.write(chunk)
    if log:
        elapsed = time.time() - t0
        log(f"  wrote {dest.name} ({dest.stat().st_size // (1 << 20)} MB, {elapsed:.1f}s)")
    return dest


def ensure_rstE3_files(*, force: bool = False, log=None) -> dict[str, Path]:
    """Download RSTE3 metadata + cellxgene expression if missing."""
    out: dict[str, Path] = {}
    for key, fname in RSTE3_FILES.items():
        dest = cache_dir() / fname
        if force and dest.is_file():
            dest.unlink()
        url = BIL_DOWNLOAD_BASE + fname
        out[key] = _download_file(url, dest, log=log)
    return out


def read_rstE3_metadata(*, force_download: bool = False, log=None) -> pd.DataFrame:
    paths = ensure_rstE3_files(force=force_download, log=log)
    df = pd.read_csv(paths["metadata"])
    if "sample_id" not in df.columns:
        raise ValueError("RSTE3 metadata missing sample_id column")
    return df


def read_rstE3_marker_expression(
    marker_genes: set[str],
    *,
    force_download: bool = False,
    log=None,
) -> pd.DataFrame:
    """
    Return genes x cells expression matrix restricted to marker_genes (case-insensitive match).
    Index = gene symbol as stored in RSTE3; columns = sample_id.
    """
    paths = ensure_rstE3_files(force=force_download, log=log)
    if log:
        log(f"  loading RSTE3 expression ({paths['cellxgene'].name})...")
    expr = pd.read_csv(paths["cellxgene"], index_col=0)
    want = {g.upper() for g in marker_genes}
    keep = [g for g in expr.index.astype(str) if g.upper() in want]
    return expr.loc[keep]
