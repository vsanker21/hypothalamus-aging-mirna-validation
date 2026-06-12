"""
Lightweight Allen Brain Cell Atlas fetch/cache (no abc_atlas_access package required).

Downloads public CSV metadata from allen-brain-cell-atlas S3 using manifest URLs.
Cache root: data/references/abc_atlas/
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

MANIFEST_URL = (
    "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/releases/20250531/manifest.json"
)

# Curated lightweight files for external validation modules.
ABC_FILES: dict[str, str] = {
    "aging_degenes": (
        "metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/aging_degenes.csv"
    ),
    "aging_cluster_mapping_pivot": (
        "metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cluster_mapping_pivot.csv"
    ),
    "aging_cell_cluster_mapping": (
        "metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cell_cluster_mapping_annotations.csv"
    ),
    "merfish_cell_metadata_cluster": (
        "metadata/MERFISH-C57BL6J-638850/20241115/views/cell_metadata_with_cluster_annotation.csv"
    ),
    "merfish_cell_metadata_parcellation": (
        "metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv"
    ),
    "merfish_ccf_coordinates": (
        "metadata/MERFISH-C57BL6J-638850-CCF/20231215/ccf_coordinates.csv"
    ),
    "merfish_gene": (
        "metadata/MERFISH-C57BL6J-638850/20241115/gene.csv"
    ),
    "merfish_example_genes": (
        "metadata/MERFISH-C57BL6J-638850/20241115/views/example_genes_all_cells_expression.csv"
    ),
}

S3_BASE = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/"


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def cache_dir() -> Path:
    return _project_root() / "data" / "references" / "abc_atlas"


def _manifest_path() -> Path:
    return cache_dir() / "manifest_20250531.json"


def load_manifest(force_refresh: bool = False) -> dict[str, Any]:
    p = _manifest_path()
    if p.is_file() and not force_refresh:
        return json.loads(p.read_text(encoding="utf-8"))
    cache_dir().mkdir(parents=True, exist_ok=True)
    r = requests.get(MANIFEST_URL, timeout=120)
    r.raise_for_status()
    data = r.json()
    p.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def download_abc_file(
    key: str,
    *,
    force: bool = False,
    log=None,
) -> Path | None:
    """Download a curated ABC file by logical key; return local path or None on failure."""
    if key not in ABC_FILES:
        raise KeyError(f"Unknown ABC file key: {key}")
    rel = ABC_FILES[key]
    dest = cache_dir() / rel.replace("/", "__")
    if dest.is_file() and not force:
        return dest
    try:
        url = S3_BASE + rel
        if log:
            log(f"ABC fetch: {key} <- {url}")
        cache_dir().mkdir(parents=True, exist_ok=True)
        for attempt in range(3):
            try:
                with requests.get(url, stream=True, timeout=300) as resp:
                    resp.raise_for_status()
                    tmp = dest.with_suffix(dest.suffix + ".part")
                    with open(tmp, "wb") as fh:
                        for chunk in resp.iter_content(chunk_size=1 << 20):
                            if chunk:
                                fh.write(chunk)
                    tmp.replace(dest)
                return dest
            except Exception as exc:
                if attempt == 2:
                    if log:
                        log(f"ABC fetch failed for {key}: {exc}")
                    return None
                time.sleep(2 ** attempt)
    except Exception as exc:
        if log:
            log(f"ABC fetch error ({key}): {exc}")
        return None
    return dest


def read_abc_csv(key: str, *, force_download: bool = False, log=None, **read_csv_kw) -> pd.DataFrame | None:
    path = download_abc_file(key, force=force_download, log=log)
    if path is None or not path.is_file():
        return None
    try:
        return pd.read_csv(path, **read_csv_kw)
    except Exception as exc:
        if log:
            log(f"ABC read failed ({key}): {exc}")
        return None
