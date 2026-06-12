"""
External HypoMap third-ventricle niche validation figure.

Requires:
  exploratory_niche_hypomap_external_validation_per_stratum.csv
  exploratory_niche_hypomap_external_validation_summary.json

Output:
  outputs/figures/sa_bundle/fig_niche_hypomap_external_validation.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def build_niche_hypomap_external_validation_figure(out_dir: Path) -> bool:
    out_dir = out_dir.resolve()
    per_path = out_dir / "exploratory_niche_hypomap_external_validation_per_stratum.csv"
    summ_path = out_dir / "exploratory_niche_hypomap_external_validation_summary.json"
    if not per_path.is_file():
        return False

    df = pd.read_csv(per_path)
    summ = json.loads(summ_path.read_text(encoding="utf-8")) if summ_path.is_file() else {}
    thr = float(summ.get("rho_threshold_hypomap_validated", 0.82))

    fig = plt.figure(figsize=(12.8, 6.0), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.32)

    # Panel A: HypoMap niche rho per stratum
    ax_a = fig.add_subplot(gs[0, 0])
    d2 = df.dropna(subset=["hypomap_niche_best_rho"]).sort_values(
        "hypomap_niche_best_rho", ascending=True
    )
    y = np.arange(len(d2))
    colors = []
    for _, row in d2.iterrows():
        if row.get("concordant_marker_and_hypomap"):
            colors.append("#1a7f4e")
        elif row.get("hypomap_niche_validated"):
            colors.append("#5dade2")
        elif row.get("is_third_ventricle_niche"):
            colors.append("#f4d03f")
        else:
            colors.append("#bdc3c7")
    ax_a.barh(y, d2["hypomap_niche_best_rho"].astype(float), color=colors, height=0.72)
    ax_a.axvline(thr, color="#c0392b", ls="--", lw=0.9, alpha=0.85, label=f"rho={thr}")
    lbl = (
        d2["stratum"].astype(str)
        + " | "
        + d2["rank1_module"].fillna("?").astype(str)
        + "\n"
        + d2["hypomap_niche_best_type"].fillna("").astype(str).str.replace("C185_", "", regex=False)
    )
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(lbl, fontsize=5.5)
    ax_a.set_xlabel("HypoMap third-ventricle niche Spearman rho (external)")
    ax_a.set_title("A  Independent atlas validation (HypoMap C185)", fontsize=10, fontweight="bold")
    ax_a.legend(loc="lower right", fontsize=7)

    # Panel B: concordance + attenuation
    ax_b = fig.add_subplot(gs[0, 1])
    groups = ["Concordant\n(marker + HypoMap)", "HypoMap only", "Marker only", "Neither"]
    conc = df["concordant_marker_and_hypomap"].astype(bool)
    hyp = df["hypomap_niche_validated"].astype(bool)
    mar = df["is_third_ventricle_niche"].astype(bool)
    masks = [
        conc,
        hyp & ~mar,
        mar & ~hyp,
        ~hyp & ~mar,
    ]
    medians = []
    counts = []
    for m in masks:
        sub = df.loc[m, "delta_median_abs_logfc"].astype(float).dropna()
        medians.append(float(sub.median()) if len(sub) else float("nan"))
        counts.append(int(len(sub)))
    xpos = np.arange(4)
    ax_b.bar(xpos, medians, color=["#1a7f4e", "#5dade2", "#f4d03f", "#bdc3c7"], width=0.55)
    ax_b.axhline(0, color="#7f8c8d", lw=0.8)
    ax_b.set_xticks(xpos)
    ax_b.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(groups, counts)], fontsize=7)
    ax_b.set_ylabel("Median Δ|logFC| (union − non-union)")
    ax_b.set_title("B  Attenuation by validation class", fontsize=10, fontweight="bold")

    n_conc = summ.get("n_concordant_strata", "?")
    fisher_p = summ.get("fisher_p_marker_vs_hypomap", float("nan"))
    jacc = summ.get("jaccard_marker_vs_hypomap_niche_ids", float("nan"))
    note = (
        f"Concordant strata: {summ.get('concordant_strata_ids', [])}\n"
        f"Jaccard(marker,HypoMap IDs)={jacc:.2f}; Fisher P={fisher_p:.3g}\n"
        f"Cluster 22 → {summ.get('cluster_22_hypomap_best_type', 'N/A')}"
    )
    ax_b.text(0.5, 0.98, note, transform=ax_b.transAxes, ha="center", va="top", fontsize=7, wrap=True)

    fig.suptitle(
        "Third-ventricle niche claim validated against independent HypoMap cell-type reference",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out_png = out_dir / "figures" / "sa_bundle" / "fig_niche_hypomap_external_validation.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    args = ap.parse_args()
    return 0 if build_niche_hypomap_external_validation_figure(args.outputs_dir) else 1


if __name__ == "__main__":
    raise SystemExit(main())
