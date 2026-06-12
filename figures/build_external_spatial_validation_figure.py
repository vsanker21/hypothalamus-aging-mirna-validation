"""Fig. 9 — External spatial/aging validation bundle (Allen Jin, MERFISH, RSTE3, ISH)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build_external_spatial_validation_figure(out_dir: Path, fig_path: Path) -> bool:
    jin = out_dir / "exploratory_allen_aging_scrna_validation_summary.json"
    mer = out_dir / "exploratory_merfish_spatial_validation_summary.json"
    rst = out_dir / "exploratory_allen_aging_spatial_validation_summary.json"
    ish = out_dir / "exploratory_allen_ish_marker_anatomy.csv"
    per_mer = out_dir / "exploratory_merfish_spatial_validation_per_stratum.csv"
    hm = out_dir / "exploratory_allen_aging_scrna_hallmark_overlap.csv"

    if not jin.is_file():
        return False

    jin_s = json.loads(jin.read_text(encoding="utf-8"))
    mer_s = json.loads(mer.read_text(encoding="utf-8")) if mer.is_file() else {}
    rst_s = json.loads(rst.read_text(encoding="utf-8")) if rst.is_file() else {}

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.5), dpi=150)

    # A — Jin V3 age-DE gene counts by cell type (top)
    de_path = out_dir / "exploratory_allen_aging_scrna_tanycyte_ependymal_agede.csv"
    ax = axes[0, 0]
    if de_path.is_file():
        de = pd.read_csv(de_path)
        counts = (
            de.groupby("grouping_name")["gene_symbol"]
            .nunique()
            .sort_values(ascending=True)
            .tail(12)
        )
        ax.barh(range(len(counts)), counts.values, color="#4a7c59")
        ax.set_yticks(range(len(counts)))
        ax.set_yticklabels([str(x)[:28] for x in counts.index], fontsize=7)
        ax.set_xlabel("# age-DE genes (Jin ABC)")
        ax.set_title("A  Jin 2025 V3 niche age-DE burden")
    else:
        ax.text(0.5, 0.5, "Jin age-DE table missing", ha="center", va="center")
        ax.set_axis_off()

    # B — Hallmark overlap with Jin V3 age-DE
    ax = axes[0, 1]
    if hm.is_file():
        h = pd.read_csv(hm).sort_values("fisher_p").head(8)
        y = np.arange(len(h))
        ax.barh(y, -np.log10(h["fisher_p"].astype(float).clip(lower=1e-300)), color="#8e44ad")
        ax.set_yticks(y)
        ax.set_yticklabels([str(t)[:32] for t in h["hallmark_term"]], fontsize=7)
        ax.set_xlabel("-log10(Fisher P) vs Jin V3 age-DE")
        ax.set_title("B  Hallmark ∩ Jin V3 aging genes")
    else:
        ax.text(0.5, 0.5, "Hallmark overlap missing", ha="center", va="center")
        ax.set_axis_off()

    # C — MERFISH spatial validation per stratum
    ax = axes[1, 0]
    if per_mer.is_file():
        pm = pd.read_csv(per_mer)
        pm = pm[pm["stratum"].notna()].copy()
        pm["stratum"] = pm["stratum"].astype(str)
        if "merfish_frac_hypothalamus" in pm.columns:
            pm = pm.sort_values("merfish_frac_hypothalamus", ascending=True)
            colors = [
                "#c0392b" if bool(v) else "#95a5a6"
                for v in pm.get("merfish_spatial_validated", pd.Series(False, index=pm.index))
            ]
            ax.barh(pm["stratum"], pm["merfish_frac_hypothalamus"].astype(float), color=colors)
            ax.set_xlabel("Fraction of MERFISH subclass cells in HY division")
            ax.set_title("C  MERFISH hypothalamus localization")
        else:
            ax.text(0.5, 0.5, "MERFISH spatial columns missing", ha="center", va="center")
            ax.set_axis_off()
    else:
        ax.text(0.5, 0.5, "MERFISH per-stratum missing", ha="center", va="center")
        ax.set_axis_off()

    # D — RSTE3 cell-level age effects or summary text
    ax = axes[1, 1]
    rst_stats_path = out_dir / "exploratory_allen_rstE3_cell_level_age_stats.csv"
    cell_level = rst_s.get("rstE3_cell_level", {})
    if rst_stats_path.is_file() and cell_level.get("spatial_data_downloaded"):
        rs = pd.read_csv(rst_stats_path)
        rs = rs.sort_values("mannwhitney_p").head(10)
        if not rs.empty:
            y = np.arange(len(rs))
            colors = ["#2980b9" if v else "#bdc3c7" for v in rs.get("direction_concordant_with_s22", [False] * len(rs))]
            ax.barh(
                y,
                rs["log2fc_aged_vs_young"].astype(float),
                color=colors,
            )
            ax.axvline(0, color="0.4", lw=0.8)
            labels = [
                f"{str(r['gene'])[:10]} ({str(r['subclass_label'])[:12]})"
                for _, r in rs.iterrows()
            ]
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("log2 FC aged vs young (RSTE3)")
            ax.set_title("D  RSTE3 V3 cell-level aging (top P)")
        else:
            ax.text(0.5, 0.5, "RSTE3 stats empty", ha="center", va="center")
            ax.set_axis_off()
    else:
        ax.set_axis_off()
        lines = [
            "External validation summary",
            f"Jin V3 age-DE genes: {jin_s.get('n_unique_v3_niche_agede_genes', 'NA')}",
            f"miRNA union ∩ Jin V3: n={jin_s.get('mirna_union_vs_jin_v3_agede_overlap', 'NA')}, "
            f"P={jin_s.get('mirna_union_vs_jin_v3_agede_fisher_p', 'NA'):.2e}"
            if isinstance(jin_s.get("mirna_union_vs_jin_v3_agede_fisher_p"), (int, float))
            else f"miRNA union ∩ Jin V3: n={jin_s.get('mirna_union_vs_jin_v3_agede_overlap', 'NA')}",
            f"MERFISH spatial-validated strata: {mer_s.get('n_strata_merfish_spatial_validated', mer_s.get('n_merfish_niche_validated_strata', 'NA'))}",
            f"RSTE3 V3 cells: {cell_level.get('n_rstE3_v3_niche_cells', 'NA')}",
            f"RSTE3 FDR<0.05 tests: {cell_level.get('n_sig_fdr05', 'NA')}",
        ]
        ax.text(0.05, 0.95, "\n".join(lines), va="top", fontsize=9, family="monospace")

    fig.suptitle("Fig. 9 — Independent spatial/aging reference validation", fontsize=11, y=1.02)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)
    return True


def build_allen_ish_marker_figure(out_dir: Path, fig_path: Path) -> bool:
    csv_path = out_dir / "exploratory_allen_ish_marker_anatomy.csv"
    if not csv_path.is_file():
        return False
    df = pd.read_csv(csv_path)
    pivot = df.pivot_table(
        index="gene", columns="structure", values="max_expression_energy", aggfunc="max"
    )
    if pivot.empty:
        return False
    fig, ax = plt.subplots(figsize=(6.5, 3.8), dpi=150)
    im = ax.imshow(pivot.astype(float).fillna(0).values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=9)
    ax.set_title("Supplementary — Allen ISH expression energy (adult anatomy)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="max energy")
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)
    return True
