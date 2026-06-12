"""Multi-library ORA (Enrichr) with optional cross-library FDR on pooled term p-values."""
from __future__ import annotations

from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# Curated for Mus musculus Enrichr (verified library names)
DEFAULT_MOUSE_ORA_LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "KEGG_2019_Mouse",
    "GO_Biological_Process_2023",
    "WikiPathways_2024_Mouse",
    "Reactome_2022",
    "Aging_Perturbations_from_GEO_up",
    "GTEx_Aging_Signatures_2021",
]


def run_multi_ora(
    genes: list[str],
    libraries: list[str] | None = None,
    organism: str = "mouse",
    top_n_terms_per_lib: int = 40,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      (per_term_long, per_library_summary)
    """
    libraries = libraries or DEFAULT_MOUSE_ORA_LIBRARIES
    rows = []
    lib_summary = []
    for lib in libraries:
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=[lib],
                organism=organism,
                outdir=None,
                no_plot=True,
            )
        except Exception as e:
            lib_summary.append({"library": lib, "status": "error", "error": str(e)})
            continue
        df = enr.results.copy()
        if len(df) == 0:
            lib_summary.append({"library": lib, "status": "empty", "n_terms": 0})
            continue
        df = df.sort_values("Adjusted P-value").head(top_n_terms_per_lib)
        df.insert(0, "library", lib)
        for _, r in df.iterrows():
            rows.append(r.to_dict())
        lib_summary.append(
            {
                "library": lib,
                "status": "ok",
                "n_terms_reported": len(df),
                "best_adj_p": float(df["Adjusted P-value"].min()),
            }
        )
    long_df = pd.DataFrame(rows)
    if len(long_df) and "P-value" in long_df.columns:
        rej, q, _, _ = multipletests(
            long_df["P-value"].astype(float),
            method="fdr_bh",
        )
        long_df["global_bh_q"] = q
        long_df["global_bh_reject_005"] = rej
    sum_df = pd.DataFrame(lib_summary)
    return long_df, sum_df


def save_multi_ora(long_df: pd.DataFrame, sum_df: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(out_dir / "enrichr_multi_library_ora_terms.csv", index=False)
    sum_df.to_csv(out_dir / "enrichr_multi_library_summary.csv", index=False)
