"""
Join Seurat-stratum cross-modal results with cluster putative labels (if present),
apply exploratory BH adjustment across strata, and write diagnostic figures.

Prerequisites (from run_extended.py + R stratified pseudobulk):
  outputs/exploratory_crossmodal_celltype_strata_summary.csv
  outputs/gse188646_strata/stratum_*_young_vs_aged_deg.csv
  outputs/exploratory_crossmodal_gene_burden_vs_aging_logfc.csv  (weighted burden per gene)
Optional:
  outputs/gse188646_cluster_annotation/cluster_putative_labels.csv

Usage (from feasibility_study/):
  python tools/diagnostic_crossmodal_strata_figures.py
  python tools/diagnostic_crossmodal_strata_figures.py --outputs-dir E:\\path\\to\\outputs
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from scipy import stats


def _bh(pvals: np.ndarray) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests

    ok = np.isfinite(pvals)
    out = np.full_like(pvals, np.nan, dtype=float)
    if ok.sum():
        _, q, _, _ = multipletests(pvals[ok], method="fdr_bh")
        out[np.where(ok)[0]] = q
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=Path(__file__).resolve().parents[1] / "outputs")
    ap.add_argument(
        "--write-supplementary-docx",
        action="store_true",
        help="After figures/CSV, run tools/generate_supplementary_figures_strata_crossmodal_docx.py",
    )
    args = ap.parse_args()
    out = args.outputs_dir.resolve()
    summ_path = out / "exploratory_crossmodal_celltype_strata_summary.csv"
    bur_path = out / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv"
    labels_path = out / "gse188646_cluster_annotation" / "cluster_putative_labels.csv"
    strata_dir = out / "gse188646_strata"
    fig_dir = out / "figures" / "crossmodal_strata_diagnostics"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not summ_path.is_file():
        print(f"Missing {summ_path.name}; run run_extended with strata DE first.", file=sys.stderr)
        return 1
    if not bur_path.is_file():
        print(f"Missing {bur_path.name}; run journal-tier cross-modal first.", file=sys.stderr)
        return 1

    d = pd.read_csv(summ_path)
    d["stratum"] = d["stratum"].astype(str).str.strip()

    if labels_path.is_file():
        lab = pd.read_csv(labels_path)
        idcol = "seurat_cluster_id" if "seurat_cluster_id" in lab.columns else lab.columns[0]
        lab[idcol] = lab[idcol].astype(str).str.strip()
        keep = [c for c in ("rank1_module", "rank2_module", "delta_z", "n_cells") if c in lab.columns]
        lab_s = lab[[idcol] + keep].rename(columns={idcol: "stratum"})
        d = d.merge(lab_s, on="stratum", how="left", suffixes=("", "_labeldup"))
    else:
        print(f"No {labels_path}; annotated table will omit putative labels.")

    for col in ("spearman_p_weighted", "perm_p_rho_gene_shuffle", "mannwhitney_abs_logFC_union_vs_nonunion_p"):
        if col in d.columns:
            d[f"{col}_fdr_bh_across_strata"] = _bh(d[col].astype(float).values)

    ann_path = out / "exploratory_crossmodal_celltype_strata_annotated.csv"
    d.to_csv(ann_path, index=False)
    print(f"Wrote {ann_path}")

    # --- Forest plot: rho per stratum ---
    d2 = d.sort_values("spearman_rho_weighted_burden_vs_logFC", ascending=True).reset_index(drop=True)
    y = np.arange(len(d2))
    rho = d2["spearman_rho_weighted_burden_vs_logFC"].astype(float).values
    if "rank1_module" in d2.columns:
        cats = pd.Categorical(d2["rank1_module"].fillna("NA"))
        colors = plt.cm.tab20(cats.codes % 20)
    else:
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(d2)))

    fig, ax = plt.subplots(figsize=(8, max(6, len(d2) * 0.18)))
    ax.barh(y, rho, color=colors, edgecolor="none", height=0.7)
    ax.axvline(0, color="0.3", lw=0.8)
    lbl = d2["stratum"].astype(str)
    if "rank1_module" in d2.columns:
        lbl = lbl + " | " + d2["rank1_module"].fillna("?").astype(str)
    ax.set_yticks(y)
    ax.set_yticklabels(lbl, fontsize=7)
    ax.set_xlabel("Spearman ρ (weighted miRNA-target burden vs stratum pseudobulk logFC)")
    ax.set_title("Cell-type / cluster–resolved cross-modal coupling (GSE188646 strata)")
    fig.tight_layout()
    fp = fig_dir / "fig_strata_spearman_rho_forest.png"
    fig.savefig(fp, dpi=160)
    plt.close(fig)
    print(f"Wrote {fp}")

    # --- Glial-focused forest (Astrocyte / Microglia in rank1 or rank2) ---
    glial_mod = frozenset({"Astrocyte", "Microglia"})

    def _glial_row(r) -> bool:
        if "rank1_module" not in d.columns:
            return False
        r1 = str(r.get("rank1_module") or "").strip()
        r2 = str(r.get("rank2_module") or "").strip()
        return r1 in glial_mod or r2 in glial_mod

    dg = d[d.apply(_glial_row, axis=1)].copy()
    fp_glial: Path | None = None
    if len(dg) >= 1:
        fp_glial = fig_dir / "fig_strata_glial_astrocyte_microglia_rho_forest.png"
        dg2 = dg.sort_values("spearman_rho_weighted_burden_vs_logFC", ascending=True).reset_index(drop=True)
        yg = np.arange(len(dg2))
        rhog = dg2["spearman_rho_weighted_burden_vs_logFC"].astype(float).values
        fig_g, axg = plt.subplots(figsize=(7, max(4, len(dg2) * 0.35)))
        col_g = np.where(
            dg2["rank1_module"].fillna("").astype(str).isin(glial_mod),
            "#2ca02c",
            np.where(dg2["rank2_module"].fillna("").astype(str).isin(glial_mod), "#9467bd", "0.45"),
        )
        axg.barh(yg, rhog, color=col_g, edgecolor="none", height=0.65)
        axg.axvline(0, color="0.3", lw=0.8)
        lblg = dg2["stratum"].astype(str) + " | " + dg2["rank1_module"].fillna("?").astype(str)
        axg.set_yticks(yg)
        axg.set_yticklabels(lblg, fontsize=8)
        axg.set_xlabel("Spearman ρ (burden vs stratum logFC)")
        axg.set_title("Strata with Astrocyte or Microglia in rank1 / rank2 (marker mapping)")
        fig_g.tight_layout()
        fig_g.savefig(fp_glial, dpi=160)
        plt.close(fig_g)
        print(f"Wrote {fp_glial}")
    else:
        print("No Astrocyte/Microglia-labelled strata; skip glial-only forest.")

    # --- Multi-panel: global pseudobulk + top |rho| strata ---
    bur = pd.read_csv(bur_path)
    bur["gene_u"] = bur["gene"].astype(str).str.strip().str.upper()
    bcols = ["gene_u", "weighted_burden"]
    if "logFC" in bur.columns:
        bcols.append("logFC")
    bur = bur[bcols]

    def panel(ax, title: str, deg: pd.DataFrame, min_pts: int = 200) -> None:
        deg = deg.copy()
        deg["gene_u"] = deg["gene"].astype(str).str.strip().str.upper()
        m = deg.merge(bur[["gene_u", "weighted_burden"]], on="gene_u", how="inner")
        if len(m) < min_pts:
            ax.text(0.5, 0.5, f"n={len(m)} (<{min_pts})", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            return
        x = m["weighted_burden"].astype(float).values
        yv = m["logFC"].astype(float).values
        ax.scatter(x, yv, s=4, alpha=0.28, c="0.25")
        r, p = stats.spearmanr(x, yv, nan_policy="omit")
        ax.set_xlabel("Weighted burden")
        ax.set_ylabel("logFC (stratum DE)")
        ax.set_title(title + f"\nρ={r:.3f}, P={p:.3g}, n={len(m)}")
        ax.axhline(0, color="0.6", lw=0.5)
        ax.axvline(0, color="0.6", lw=0.5)

    global_deg = out / "gse188646_young_vs_aged_deg.csv"
    if not global_deg.is_file():
        print("Missing gse188646_young_vs_aged_deg.csv; skip scatter multipanel.", file=sys.stderr)
        return 0

    top = d.assign(_ar=d["spearman_rho_weighted_burden_vs_logFC"].abs()).sort_values("_ar", ascending=False).head(3)

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.28)
    ax0 = fig.add_subplot(gs[0, 0])
    gdeg = pd.read_csv(global_deg)
    panel(ax0, "Global pseudobulk (all nuclei)", gdeg)
    positions = [(0, 1), (1, 0), (1, 1)]
    for ax_idx, (_, row) in enumerate(top.iterrows()):
        if ax_idx >= 3:
            break
        r, c = positions[ax_idx]
        ax = fig.add_subplot(gs[r, c])
        st = str(row["stratum"]).strip()
        fn = row.get("deg_csv") or f"stratum_{st}_young_vs_aged_deg.csv"
        pth = strata_dir / str(fn)
        if not pth.is_file():
            ax.text(0.5, 0.5, f"missing {fn}", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"Stratum {st}")
            continue
        sdeg = pd.read_csv(pth)
        mod = ""
        if "rank1_module" in row and pd.notna(row["rank1_module"]):
            mod = f" ({row['rank1_module']})"
        panel(ax, f"Stratum {st}{mod}", sdeg)

    fig.suptitle("Burden vs logFC: global + top 3 strata by |Spearman ρ|", fontsize=12, y=1.02)
    fp2 = fig_dir / "fig_strata_burden_logfc_scatter_grid.png"
    fig.savefig(fp2, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {fp2}")

    figs = [
        "figures/crossmodal_strata_diagnostics/fig_strata_spearman_rho_forest.png",
        "figures/crossmodal_strata_diagnostics/fig_strata_burden_logfc_scatter_grid.png",
    ]
    if fp_glial is not None and fp_glial.is_file():
        figs.append("figures/crossmodal_strata_diagnostics/fig_strata_glial_astrocyte_microglia_rho_forest.png")

    mini = {
        "n_strata": int(len(d)),
        "annotated_csv": "exploratory_crossmodal_celltype_strata_annotated.csv",
        "figures": figs,
        "note": "FDR_BH across strata is exploratory multiplicity control on independent cluster tests.",
    }
    (fig_dir / "README_DIAGNOSTICS.json").write_text(json.dumps(mini, indent=2), encoding="utf-8")

    if args.write_supplementary_docx:
        proj = Path(__file__).resolve().parents[1]
        gen = proj / "tools" / "generate_supplementary_figures_strata_crossmodal_docx.py"
        if gen.is_file():
            r = subprocess.run(
                [sys.executable, str(gen), "--outputs-dir", str(out)],
                cwd=str(proj),
                capture_output=True,
                text=True,
                timeout=120,
            )
            if r.stdout:
                print(r.stdout.strip())
            if r.returncode != 0:
                print(r.stderr or "", file=sys.stderr)
                return int(r.returncode)
        else:
            print(f"Missing {gen}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
