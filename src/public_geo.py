"""GEO metadata for public hypothalamus aging datasets (no large matrix download)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

try:
    from GEOparse import get_GEO
except ImportError:
    get_GEO = None

from config import GEO_HYPOTHALAMUS_AGING, OUTPUT_DIR, PROJECT_ROOT


def fetch_gse_sample_table(gse_id: str, destdir: Path | None = None) -> pd.DataFrame:
    if get_GEO is None:
        raise ImportError("Install GEOparse: pip install GEOparse")
    destdir = destdir or (PROJECT_ROOT / "data" / "geo_cache")
    destdir.mkdir(parents=True, exist_ok=True)
    gse = get_GEO(geo=gse_id, destdir=str(destdir), silent=True)
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        ch = gsm.metadata.get("characteristics_ch1", [])
        title = (gsm.metadata.get("title") or [""])[0]
        src = (gsm.metadata.get("source_name_ch1") or [""])[0]
        rows.append(
            {
                "gsm": gsm_name,
                "title": title,
                "source_name_ch1": src,
                "characteristics_ch1": " | ".join(ch) if isinstance(ch, list) else str(ch),
            }
        )
    return pd.DataFrame(rows)


def save_geo_context(gse_id: str = GEO_HYPOTHALAMUS_AGING) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = fetch_gse_sample_table(gse_id)
    out = OUTPUT_DIR / f"{gse_id}_sample_metadata.csv"
    df.to_csv(out, index=False)
    meta = {
        "gse": gse_id,
        "description": (
            "Single-nucleus RNA-seq of hypothalamic nuclei from young vs aged female mice; "
            "40,064 nuclei in original publication. Use for contextual rationale: hypothalamic "
            "aging is resolvable at single-cell resolution; complements bulk miRNA/exosome "
            "mechanisms in Zhang et al. (2017)."
        ),
        "reference": "PMID 34782725 / Nat Aging 2022 (dataset GSE188646)",
        "note": (
            "Raw integrated object is distributed as a large RDS file on GEO; this pipeline "
            "records sample metadata only to avoid multi-GB downloads while preserving "
            "citable public-data linkage."
        ),
    }
    with open(OUTPUT_DIR / f"{gse_id}_context.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    return out
