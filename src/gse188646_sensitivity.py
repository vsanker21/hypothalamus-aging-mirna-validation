"""Exploratory GSE188646 DE overlap sensitivity (Fisher) + DE calibration table."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import overlap_fisher


def run_gse188646_fisher_sensitivity(
    union: set[str],
    uni_mmu: set[str],
    gse_deg: Path,
    out_dir: Path,
    log,
) -> None:
    """
    Primary Fisher in run_extended uses padj<=0.05. This block documents calibration and
    runs labelled exploratory overlaps (padj 0.1; top-200 / top-500 rank-score gene sets).
    """
    if not gse_deg.is_file():
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    cal = overlap_fisher.gse188646_de_calibration(gse_deg)
    cal.to_csv(out_dir / "exploratory_gse188646_de_padj_calibration.csv", index=False)
    log("\n=== Exploratory: GSE188646 pseudobulk DE calibration (padj thresholds) ===")
    log(cal.to_string(index=False))

    uni_u = {g.upper() for g in uni_mmu}

    def _fisher(label: str, ext: set[str], suffix: str) -> None:
        ext_in_u = {g.upper() for g in ext} & uni_u
        res = overlap_fisher.fisher_overlap(union, ext_in_u, uni_mmu)
        row = {"label": label, **res}
        pd.DataFrame([row]).to_csv(out_dir / f"fisher_targets_vs_gse188646_de_{suffix}.csv", index=False)
        log(f"Exploratory Fisher [{label}]: {res}")

    # padj 0.1
    ext01 = overlap_fisher._read_seurat_markers(gse_deg, padj_max=0.1)
    _fisher("padj<=0.1", ext01, "exploratory_padj0.1")

    for n in (200, 500):
        ext_top = overlap_fisher.deg_top_by_rank_score(gse_deg, top_n=n)
        _fisher(f"top{n}_by_signlogFC_neglog10p", ext_top, f"exploratory_top{n}_ranked")

    note = out_dir / "exploratory_GSE188646_FISHER_SENSITIVITY_README.txt"
    note.write_text(
        "EXPLORATORY ONLY — not confirmatory.\n"
        "- Primary overlap test in EXTENDED_REPORT uses padj<=0.05 on pseudobulk edgeR output.\n"
        "- This folder adds padj<=0.1 and rank-based top gene sets for sensitivity when FDR is sparse.\n"
        "- Do not upgrade these to primary claims without pre-registration and independent replication.\n",
        encoding="utf-8",
    )
    log(f"Wrote {note.name}")
