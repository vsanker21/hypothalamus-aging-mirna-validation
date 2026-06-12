"""
Publication figure: third-ventricle niche localization of miRNA-target aging lability.

Requires:
  exploratory_niche_lability_per_stratum.csv
  exploratory_niche_lability_localization_summary.json

Output:
  outputs/figures/sa_bundle/fig_niche_lability_localization.png
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def build_niche_lability_figure(out_dir: Path) -> bool:
    out_dir = out_dir.resolve()
    per_path = out_dir / "exploratory_niche_lability_per_stratum.csv"
    summ_path = out_dir / "exploratory_niche_lability_localization_summary.json"
    if not per_path.is_file():
        return False

    df = pd.read_csv(per_path)
    summ = json.loads(summ_path.read_text(encoding="utf-8")) if summ_path.is_file() else {}

    niche_mask = df["is_third_ventricle_niche"].astype(bool)
    colors = np.where(niche_mask, "#1a7f4e", "#7f8c8d")

    fig = plt.figure(figsize=(12.5, 6.2), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.32)

    # Panel A: forest of -log10(MW p) for lability, sorted
    ax_a = fig.add_subplot(gs[0, 0])
    d2 = df.sort_values("mannwhitney_abs_logFC_union_vs_nonunion_p", ascending=True).reset_index(drop=True)
    y = np.arange(len(d2))
    pvals = d2["mannwhitney_abs_logFC_union_vs_nonunion_p"].astype(float).clip(lower=1e-300)
    x = -np.log10(pvals)
    ncol = np.where(d2["is_third_ventricle_niche"].astype(bool), "#1a7f4e", "#7f8c8d")
    ax_a.barh(y, x, color=ncol, height=0.72, edgecolor="none")
    ax_a.axvline(-np.log10(0.05), color="#c0392b", ls="--", lw=0.9, alpha=0.8)
    lbl = (
        d2["stratum"].astype(str)
        + " | "
        + d2["rank1_module"].fillna("?").astype(str)
        + ("+" + d2["rank2_module"].fillna("").astype(str)).where(d2["rank2_module"].notna(), "")
    )
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(lbl, fontsize=6.5)
    ax_a.set_xlabel(r"$-\log_{10}$(Mann–Whitney P) for |logFC| union vs non-union")
    ax_a.set_title("A  Per-cluster lability (cell-type–resolved pseudobulk)", fontsize=10, fontweight="bold")

    # Panel B: delta median |logFC| niche vs other + annotation
    ax_b = fig.add_subplot(gs[0, 1])
    groups = ["Third-ventricle\nniche panel", "Other\nclusters"]
    med_n = summ.get("median_delta_abs_logfc_niche")
    med_o = summ.get("median_delta_abs_logfc_other")
    wm_n = summ.get("weighted_mean_delta_abs_logfc_niche")
    wm_o = summ.get("weighted_mean_delta_abs_logfc_other")
    if med_n is None or med_o is None:
        med_n = float(df.loc[niche_mask, "delta_median_abs_logfc"].median()) if niche_mask.any() else 0.0
        med_o = float(df.loc[~niche_mask, "delta_median_abs_logfc"].median()) if (~niche_mask).any() else 0.0
    ax_b.bar([0, 1], [med_n, med_o], color=["#1a7f4e", "#7f8c8d"], width=0.55, alpha=0.9)
    ax_b.set_xticks([0, 1])
    ax_b.set_xticklabels(groups)
    ax_b.set_ylabel("Median Δ|logFC| (union − non-union)")
    ax_b.set_title("B  Niche vs non-niche strata (median effect size)", fontsize=10, fontweight="bold")
    wp = summ.get("wilcoxon_p_delta_median_niche_vs_other", float("nan"))
    pp = summ.get("perm_p_median_delta_diff_ge_obs", float("nan"))
    gp = summ.get("global_pseudobulk_mannwhitney_p_reference")
    note = (
        f"Niche strata n={summ.get('n_third_ventricle_niche_strata', niche_mask.sum())}; "
        f"Wilcoxon P={wp:.3g}; perm P={pp:.3g}\n"
        f"Weighted mean Δ|logFC|: niche={wm_n:.4f}, other={wm_o:.4f}\n"
        f"Global pseudobulk MW P (ref)={gp}"
    )
    ax_b.text(0.5, 0.98, note, transform=ax_b.transAxes, ha="center", va="top", fontsize=8, wrap=True)

    fig.suptitle(
        "MiRNA-target aging lability is least attenuated in third-ventricle tanycyte/NSC-niche clusters (GSE188646 strata)",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out_png = out_dir / "figures" / "sa_bundle" / "fig_niche_lability_localization.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    args = ap.parse_args()
    return 0 if build_niche_lability_figure(args.outputs_dir) else 1


if __name__ == "__main__":
    raise SystemExit(main())
