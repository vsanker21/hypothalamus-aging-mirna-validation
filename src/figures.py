"""Figures for feasibility outputs."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from config import OUTPUT_DIR


def setup_style():
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["figure.dpi"] = 120


def plot_mirna_volcano(df: Path | pd.DataFrame, out: Path):
    if isinstance(df, Path):
        d = pd.read_csv(df)
    else:
        d = df.copy()
    setup_style()
    fig, ax = plt.subplots(figsize=(8, 6))
    x = d["logfc_htnsc_vs_astro"]
    y = np.abs(x)
    ax.scatter(x, y, alpha=0.35, s=12, c="#2c7fb8", edgecolors="none")
    top = d.nlargest(25, "logfc_htnsc_vs_astro")
    ax.scatter(top["logfc_htnsc_vs_astro"], np.abs(top["logfc_htnsc_vs_astro"]), c="#e34a33", s=40)
    for _, r in top.head(8).iterrows():
        ax.text(
            r["logfc_htnsc_vs_astro"],
            abs(r["logfc_htnsc_vs_astro"]),
            str(r["mirna"]).split("-")[-1][:10],
            fontsize=7,
            ha="left",
        )
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_xlabel("log expression difference (htNSC mean − astrocyte mean)")
    ax.set_ylabel("|Δ| (ranking emphasis)")
    ax.set_title("miRNA enrichment toward htNSCs vs astrocytes (MOESM2)")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_cytokine_paired(df: pd.DataFrame, out: Path):
    setup_style()
    rows = []
    for _, r in df.iterrows():
        rows.append({"gene": r["gene"], "arm": "Vehicle", "value": r["vehicle"]})
        rows.append({"gene": r["gene"], "arm": "Exosome", "value": r["exosome"]})
    pd_long = pd.DataFrame(rows)
    g = sns.catplot(
        data=pd_long,
        x="arm",
        y="value",
        hue="arm",
        col="gene",
        kind="box",
        sharey=False,
        height=3,
        aspect=0.9,
        palette="Set2",
        legend=False,
    )
    g.set_axis_labels("", "mRNA (normalized)")
    g.fig.subplots_adjust(top=0.85)
    g.fig.suptitle("Exosome vs vehicle (paired animals, MOESM22)")
    out.parent.mkdir(parents=True, exist_ok=True)
    g.savefig(out, bbox_inches="tight")
    plt.close(g.fig)


def plot_sox2_paired(df: pd.DataFrame, out: Path):
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    for ax, (a, b, title) in zip(
        axes,
        [
            ("sox2_3v_con", "sox2_3v_tk", "Sox2+ cells — 3V wall"),
            ("sox2_arc_con", "sox2_arc_tk", "Sox2+ cells — ARC"),
        ],
    ):
        x0 = np.zeros(len(df))
        x1 = np.ones(len(df))
        for i, row in df.iterrows():
            ax.plot([0, 1], [row[a], row[b]], color="#636363", alpha=0.6, lw=1)
        ax.scatter(x0, df[a], color="#3182bd", s=45, zorder=3)
        ax.scatter(x1, df[b], color="#de2d26", s=45, zorder=3)
        ax.set_xticks([0, 1], ["Control", "TK/GCV"])
        ax.set_title(title)
    axes[0].set_ylabel("Cell counts")
    fig.suptitle("Hypothalamic NSC-associated Sox2+ population loss with ablation (MOESM10)")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
