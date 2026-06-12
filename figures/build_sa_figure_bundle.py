"""
Build publication-style PNG panels from outputs listed in data/provenance/FIGURE_PANEL_MANIFEST.csv.

Run from feasibility_study/:
  python figures/build_sa_figure_bundle.py
  python figures/build_sa_figure_bundle.py --outputs-dir E:\\path\\to\\outputs

Optional at end of run_extended.py unless SKIP_SA_FIGURE_BUNDLE=1.
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
DEFAULT_MANIFEST = PROJECT_ROOT / "data" / "provenance" / "FIGURE_PANEL_MANIFEST.csv"


def _safe_read_csv(p: Path) -> pd.DataFrame | None:
    if not p.is_file():
        return None
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def plot_scatter_burden_logfc(csv_path: Path, title: str, out_path: Path) -> bool:
    df = _safe_read_csv(csv_path)
    if df is None or "weighted_burden" not in df.columns or "logFC" not in df.columns:
        return False
    x = df["weighted_burden"].astype(float)
    y = df["logFC"].astype(float)
    from scipy import stats

    r, p = stats.spearmanr(x, y, nan_policy="omit")
    rho = float(getattr(r, "statistic", getattr(r, "correlation", np.nan)))
    pv = float(getattr(r, "pvalue", np.nan))
    fig, ax = plt.subplots(figsize=(5.2, 4.6), dpi=150)
    ax.scatter(x, y, s=8, alpha=0.35, c="#1f4e79", edgecolors="none")
    ax.set_xlabel("htNSC-weighted miRNA target burden")
    ax.set_ylabel("GSE188646 pseudobulk logFC (Aged vs Young)")
    ax.set_title(title)
    ax.text(
        0.02,
        0.98,
        f"Spearman ρ = {rho:.3f}\np = {pv:.2e}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )
    ax.axhline(0, color="gray", lw=0.6, ls="--")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_hist_perm_rho(
    perm_csv: Path, summary_json: Path | None, title: str, out_path: Path
) -> bool:
    perm = _safe_read_csv(perm_csv)
    if perm is None or perm.empty:
        return False
    col = perm.columns[0]
    nulls = perm[col].astype(float).values
    obs = np.nan
    if summary_json is not None and summary_json.is_file():
        try:
            summ = json.loads(summary_json.read_text(encoding="utf-8"))
            obs = float(summ.get("spearman_rho_weighted_burden_vs_logFC", np.nan))
        except Exception:
            pass
    fig, ax = plt.subplots(figsize=(5.2, 4.2), dpi=150)
    ax.hist(nulls, bins=40, color="#6b8e9f", edgecolor="white", alpha=0.9)
    if np.isfinite(obs):
        ax.axvline(obs, color="#c0392b", lw=2, label=f"Observed ρ = {obs:.3f}")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel("Permutation Spearman ρ (gene-label shuffle)")
    ax.set_ylabel("Count")
    ax.set_title(title)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_volcano(deg_csv: Path, title: str, out_path: Path) -> bool:
    df = _safe_read_csv(deg_csv)
    if df is None or "logFC" not in df.columns:
        return False
    padj_col = "p_val_adj" if "p_val_adj" in df.columns else ("padj" if "padj" in df.columns else None)
    if padj_col is None:
        return False
    x = df["logFC"].astype(float)
    y = -np.log10(np.maximum(df[padj_col].astype(float), 1e-300))
    sig = df[padj_col].astype(float) < 0.05
    fig, ax = plt.subplots(figsize=(5.2, 4.6), dpi=150)
    ax.scatter(x[~sig], y[~sig], s=6, alpha=0.25, c="gray", edgecolors="none")
    ax.scatter(x[sig], y[sig], s=8, alpha=0.45, c="#1f4e79", edgecolors="none")
    ax.set_xlabel("logFC (Aged vs Young)")
    ax.set_ylabel(f"-log10({padj_col})")
    ax.set_title(title)
    ax.axvline(0, color="gray", lw=0.5)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_scatter_burden_beta_dl(csv_path: Path, title: str, out_path: Path) -> bool:
    df = _safe_read_csv(csv_path)
    if df is None or "weighted_burden" not in df.columns or "beta_DL" not in df.columns:
        return False
    x = df["weighted_burden"].astype(float)
    y = df["beta_DL"].astype(float)
    from scipy import stats

    r, p = stats.spearmanr(x, y, nan_policy="omit")
    rho = float(getattr(r, "statistic", getattr(r, "correlation", np.nan)))
    pv = float(getattr(r, "pvalue", np.nan))
    fig, ax = plt.subplots(figsize=(5.2, 4.6), dpi=150)
    ax.scatter(x, y, s=8, alpha=0.35, c="#6a3d9a", edgecolors="none")
    ax.set_xlabel("htNSC-weighted miRNA target burden")
    ax.set_ylabel("Two-cohort DL meta β (exploratory)")
    ax.set_title(title)
    ax.text(
        0.02,
        0.98,
        f"Spearman ρ = {rho:.3f}\np = {pv:.2e}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )
    ax.axhline(0, color="gray", lw=0.6, ls="--")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_violin_abs_logfc_union(
    deg_csv: Path,
    targets_long_csv: Path,
    title: str,
    out_path: Path,
    *,
    summary_json: Path | None = None,
) -> bool:
    """Mann–Whitney |logFC|: miRNA-target union vs non-union on full DE table."""
    deg = _safe_read_csv(deg_csv)
    long_df = _safe_read_csv(targets_long_csv)
    if deg is None or "logFC" not in deg.columns or long_df is None:
        return False
    tcol = "target_gene" if "target_gene" in long_df.columns else None
    if tcol is None:
        return False
    union = {str(g).strip().upper() for g in long_df[tcol].astype(str)}
    deg = deg.copy()
    deg["gene_u"] = deg["gene"].astype(str).str.strip().str.upper() if "gene" in deg.columns else deg.iloc[:, 0].astype(str).str.upper()
    in_u = deg["gene_u"].isin(union)
    u_fc = np.abs(deg.loc[in_u, "logFC"].astype(float).values)
    o_fc = np.abs(deg.loc[~in_u, "logFC"].astype(float).values)
    if len(u_fc) < 10 or len(o_fc) < 10:
        return False
    from scipy import stats

    mw = stats.mannwhitneyu(u_fc, o_fc, alternative="two-sided")
    mw_p = float(mw.pvalue)
    med_u, med_o = float(np.median(u_fc)), float(np.median(o_fc))
    fig, ax = plt.subplots(figsize=(5.0, 4.6), dpi=150)
    parts = ax.violinplot(
        [u_fc, o_fc],
        positions=[1, 2],
        showmeans=True,
        showmedians=True,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#1f4e79")
        pc.set_alpha(0.35)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["miRNA-target\nunion", "Non-union\n(DE table)"])
    ax.set_ylabel("|logFC| (Aged vs Young)")
    ax.set_title(title)
    less_p = float("nan")
    if summary_json is not None and summary_json.is_file():
        try:
            summ = json.loads(summary_json.read_text(encoding="utf-8"))
            rho = float(summ.get("spearman_rho_weighted_burden_vs_logFC", np.nan))
            rp = float(summ.get("spearman_p_weighted", np.nan))
            less_p = float(summ.get("mannwhitney_abs_logFC_union_less_p", np.nan))
            if np.isfinite(rho):
                rho_txt = f"\nSpearman(burden, logFC) ρ={rho:.3f}, p={rp:.2e}"
        except Exception:
            pass
    less_txt = f"\nMW one-sided (union < non-union) P={less_p:.2e}" if np.isfinite(less_p) else ""
    ax.text(
        0.02,
        0.98,
        f"Mann–Whitney (two-sided) P = {mw_p:.2e}\nmedian |logFC| union={med_u:.3f}\nnon-union={med_o:.3f}{less_txt}{rho_txt}",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_bar_sa_nsc_nes(csv_path: Path, title: str, out_path: Path) -> bool:
    df = _safe_read_csv(csv_path)
    if df is None or "NES" not in df.columns:
        return False
    term_col = "Term" if "Term" in df.columns else df.columns[0]
    terms = df[term_col].astype(str)
    nes = df["NES"].astype(float)
    order = np.argsort(nes.values)
    fig, ax = plt.subplots(figsize=(6.2, 4.8), dpi=150)
    colors = np.where(nes.values[order] >= 0, "#c0392b", "#2980b9")
    ax.barh(terms.values[order], nes.values[order], color=colors, alpha=0.85)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("NES (prerank GSEA vs curated niche priors)")
    ax.set_title(title)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description="Build SA-oriented figure bundle from manifest.")
    ap.add_argument(
        "--outputs-dir",
        type=Path,
        default=PROJECT_ROOT / "outputs",
        help="Directory containing pipeline CSV/JSON outputs",
    )
    ap.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="CSV manifest (panel_id, plot_type, title, ...)",
    )
    args = ap.parse_args()
    out_root: Path = args.outputs_dir.resolve()
    man_path: Path = args.manifest.resolve()
    if not man_path.is_file():
        print(f"Manifest not found: {man_path}", file=sys.stderr)
        return 1
    man = pd.read_csv(man_path)
    required = {"panel_id", "plot_type", "title", "input_csv", "out_rel"}
    if not required.issubset(set(man.columns)):
        print(f"Manifest missing columns; need {required}", file=sys.stderr)
        return 1
    ok = 0
    skip = 0
    for _, row in man.iterrows():
        plot_type = str(row["plot_type"]).strip()
        title = str(row["title"]).strip()
        out_rel = str(row["out_rel"]).strip()
        out_path = (out_root / out_rel).resolve()
        in1_raw = str(row["input_csv"]).strip() if pd.notna(row["input_csv"]) else ""
        in1 = (out_root / in1_raw) if in1_raw else None
        in2_raw = row["input_csv_secondary"] if "input_csv_secondary" in man.columns else None
        if in2_raw is not None and pd.notna(in2_raw) and str(in2_raw).strip():
            in2p: Path | None = out_root / str(in2_raw).strip()
        else:
            in2p = None

        built = False
        if plot_type == "scatter_burden_logfc":
            built = in1 is not None and plot_scatter_burden_logfc(in1, title, out_path)
        elif plot_type == "hist_perm_rho":
            built = in1 is not None and plot_hist_perm_rho(in1, in2p, title, out_path)
        elif plot_type == "volcano_logfc_padj":
            built = in1 is not None and plot_volcano(in1, title, out_path)
        elif plot_type == "scatter_burden_beta_dl":
            built = in1 is not None and plot_scatter_burden_beta_dl(in1, title, out_path)
        elif plot_type == "bar_sa_nsc_nes":
            built = in1 is not None and plot_bar_sa_nsc_nes(in1, title, out_path)
        elif plot_type == "violin_abs_logfc_union":
            summ = out_root / "exploratory_crossmodal_mirna_aging_summary.json"
            built = (
                in1 is not None
                and in2p is not None
                and plot_violin_abs_logfc_union(in1, in2p, title, out_path, summary_json=summ)
            )
        elif plot_type == "pathway_convergence":
            import importlib.util

            pc_path = PROJECT_ROOT / "figures" / "build_pathway_convergence_figure.py"
            spec = importlib.util.spec_from_file_location("build_pathway_convergence_figure", pc_path)
            if spec is None or spec.loader is None:
                built = False
            else:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                built = bool(mod.build_pathway_convergence_figure(out_root))
        elif plot_type == "niche_lability_localization":
            import importlib.util

            nl_path = PROJECT_ROOT / "figures" / "build_niche_lability_figure.py"
            spec = importlib.util.spec_from_file_location("build_niche_lability_figure", nl_path)
            if spec is None or spec.loader is None:
                built = False
            else:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                built = bool(mod.build_niche_lability_figure(out_root))
        elif plot_type == "niche_hypomap_external_validation":
            import importlib.util

            nh_path = PROJECT_ROOT / "figures" / "build_niche_hypomap_external_validation_figure.py"
            spec = importlib.util.spec_from_file_location(
                "build_niche_hypomap_external_validation_figure", nh_path
            )
            if spec is None or spec.loader is None:
                built = False
            else:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                built = bool(mod.build_niche_hypomap_external_validation_figure(out_root))
        elif plot_type == "crosscohort_lability_replication":
            import importlib.util

            cc_path = PROJECT_ROOT / "figures" / "build_crosscohort_lability_figure.py"
            spec = importlib.util.spec_from_file_location("build_crosscohort_lability_figure", cc_path)
            if spec is None or spec.loader is None:
                built = False
            else:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                built = bool(mod.build_crosscohort_lability_figure(out_root))
        elif plot_type == "external_spatial_validation":
            import importlib.util

            es_path = PROJECT_ROOT / "figures" / "build_external_spatial_validation_figure.py"
            spec = importlib.util.spec_from_file_location(
                "build_external_spatial_validation_figure", es_path
            )
            if spec is None or spec.loader is None:
                built = False
            else:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                built = bool(
                    mod.build_external_spatial_validation_figure(
                        out_root, out_root / "figures" / "sa_bundle" / "fig_external_spatial_validation.png"
                    )
                )
        elif plot_type == "allen_ish_marker_anatomy":
            import importlib.util

            es_path = PROJECT_ROOT / "figures" / "build_external_spatial_validation_figure.py"
            spec = importlib.util.spec_from_file_location(
                "build_external_spatial_validation_figure", es_path
            )
            if spec is None or spec.loader is None:
                built = False
            else:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                built = bool(
                    mod.build_allen_ish_marker_figure(
                        out_root, out_root / "figures" / "sa_bundle" / "fig_allen_ish_marker_anatomy.png"
                    )
                )
        else:
            print(f"Unknown plot_type {plot_type} for {row.get('panel_id')}", file=sys.stderr)
            skip += 1
            continue
        if built:
            ok += 1
            print(f"OK  {row['panel_id']} -> {out_path}")
        else:
            skip += 1
            print(f"SKIP {row['panel_id']} (missing data or columns)")
    print(f"Figure bundle: {ok} written, {skip} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
