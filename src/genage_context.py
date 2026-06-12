"""Download GenAge model-organism gene list (HAGR) and overlap with miRNA targets."""
from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pandas as pd
import requests

from gmt_io import CACHE_DIR

GENAGE_MODELS_ZIP = "https://genomics.senescence.info/genes/models_genes.zip"
GENAGE_CACHE_ZIP = CACHE_DIR / "genage_models_genes.zip"


def ensure_genage_models_table() -> pd.DataFrame:
    if not GENAGE_CACHE_ZIP.exists():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with requests.get(GENAGE_MODELS_ZIP, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(GENAGE_CACHE_ZIP, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    if chunk:
                        f.write(chunk)
    zf = zipfile.ZipFile(GENAGE_CACHE_ZIP)
    with zf.open("genage_models.csv") as fh:
        return pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8", errors="replace"))


def mouse_genage_symbols() -> set[str]:
    df = ensure_genage_models_table()
    if "organism" not in df.columns or "symbol" not in df.columns:
        raise ValueError("Unexpected GenAge CSV schema")
    m = df[df["organism"].astype(str).str.contains("Mus musculus", case=False, na=False)]
    syms = {str(s).strip() for s in m["symbol"].dropna() if str(s).strip()}
    return syms
