"""
Cross-cohort replication figure for aging lability (Mann–Whitney |effect|).

Requires exploratory_crosscohort_lability_replication.csv

Output: outputs/figures/sa_bundle/fig_crosscohort_lability_replication.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]

_LABELS = {
    "GSE188646_pseudobulk": "GSE188646\n(female snRNA pseudobulk)",
    "GSE87102_microarray": "GSE87102\n(male microarray)",
    "two_cohort_DL_meta": "Two-cohort\nDL meta |beta|",
}


def build_crosscohort_lability_figure(out_dir: Path) -> bool:
    out_dir = out_dir.resolve()
    csv_path = out_dir / "exploratory_crosscohort_lability_replication.csv"
    if not csv_path.is_file():
        return False

    df = pd.read_csv(csv_path)
    order = ["GSE188646_pseudobulk", "GSE87102_microarray", "two_cohort_DL_meta"]
    df["cohort_label"] = pd.Categorical(df["cohort_label"], categories=order, ordered=True)
    df = df.sort_values("cohort_label").reset_index(drop=True)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), dpi=150)

    # Panel A: -log10(MW P)
    ax = axes[0]
    pvals = df["mannwhitney_abs_effect_union_vs_nonunion_p"].astype(float).clip(lower=1e-300)
    x = -np.log10(pvals)
    colors = ["#2c3e50", "#1a5276", "#117a65"][: len(df)]
    bars = ax.bar(range(len(df)), x, color=colors, width=0.55)
    ax.axhline(-np.log10(0.05), color="#c0392b", ls="--", lw=0.9, alpha=0.85)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels([_LABELS.get(str(c), str(c)) for c in df["cohort_label"]], fontsize=8)
    ax.set_ylabel(r"$-\log_{10}$(Mann–Whitney P)")
    ax.set_title("A  Two-sided Mann–Whitney P by cohort", fontsize=10, fontweight="bold")
    for i, (bar, p) in enumerate(zip(bars, pvals)):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.15,
            f"P={p:.2e}" if p < 0.001 else f"P={p:.3g}",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    # Panel B: delta median |effect|
    ax2 = axes[1]
    deltas = df["delta_median_abs_effect"].astype(float)
    ax2.bar(range(len(df)), deltas, color=colors, width=0.55)
    ax2.axhline(0, color="#7f8c8d", lw=0.8)
    ax2.set_xticks(range(len(df)))
    ax2.set_xticklabels([_LABELS.get(str(c), str(c)) for c in df["cohort_label"]], fontsize=8)
    ax2.set_ylabel("Median Δ|effect| (union − non-union)")
    ax2.set_title("B  Median Δ|effect| (union − non-union)", fontsize=10, fontweight="bold")

    fig.suptitle(
        "Cross-cohort sensitivity: |effect| distribution in miRNA-target union vs non-union",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out_png = out_dir / "figures" / "sa_bundle" / "fig_crosscohort_lability_replication.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    args = ap.parse_args()
    return 0 if build_crosscohort_lability_figure(args.outputs_dir) else 1


if __name__ == "__main__":
    raise SystemExit(main())
