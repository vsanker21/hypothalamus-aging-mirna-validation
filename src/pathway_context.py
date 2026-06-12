"""Optional pathway context via g:Profiler (REST) — gene-symbol lists only."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

from config import GPROFILER_URL, OUTPUT_DIR


def gprofiler_ora(
    genes: Iterable[str],
    sources: tuple[str, ...] = ("GO:BP", "REAC", "KEGG"),
    organism: str = "mmusculus",
) -> pd.DataFrame:
    query = [g.strip() for g in genes if g and str(g).strip()]
    payload = {
        "query": query,
        "organism": organism,
        "sources": list(sources),
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "all_results": False,
        "background": [],
    }
    r = requests.post(GPROFILER_URL, json=payload, timeout=120)
    r.raise_for_status()
    data = r.json()
    if not data.get("result"):
        return pd.DataFrame()
    rows = []
    for row in data["result"]:
        rows.append(
            {
                "source": row.get("source"),
                "native": row.get("native"),
                "name": row.get("name"),
                "p_value": row.get("p_value"),
                "term_size": row.get("term_size"),
                "query_size": row.get("query_size"),
                "intersection_size": row.get("intersection_size"),
            }
        )
    return pd.DataFrame(rows)


def save_pathway_results(df: pd.DataFrame, fname: str) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / fname
    df.to_csv(path, index=False)
    return path


def inflammation_and_nsc_query_genes() -> list[str]:
    """
    Curated mouse gene symbols tying Zhang exosome readouts to NSC / inflammatory context.
    Small query set — suitable for ORA interpretation, not a full transcriptome substitute.
    """
    return [
        "Sox2",
        "Bmi1",
        "Nestin",
        "Msi1",
        "Cxcr4",
        "Nfkb1",
        "Rela",
        "Ikbkb",
        "Mapk1",
        "Mtor",
        "Il6",
        "Tnf",
        "Il1b",
        "Cxcl10",
        "Gfap",
        "Aif1",
        "Dcx",
        "Ascl1",
    ]
