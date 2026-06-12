"""Additional public GEO series: sample metadata tables (no large matrix downloads)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from config import OUTPUT_DIR, PROJECT_ROOT
from public_geo import fetch_gse_sample_table

# Hypothalamus / brain-related public series for contextual linkage (not merged statistically with MOESM)
EXTRA_GSE_METADATA_IDS = [
    "GSE94940",
    "GSE161219",
]


def save_extra_geo_sample_tables(out_dir: Path | None = None) -> list[Path]:
    out_dir = out_dir or OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for gse in EXTRA_GSE_METADATA_IDS:
        try:
            df = fetch_gse_sample_table(gse, destdir=PROJECT_ROOT / "data" / "geo_cache")
            p = out_dir / f"{gse}_extra_sample_metadata.csv"
            df.to_csv(p, index=False)
            written.append(p)
        except Exception:
            continue
    return written
