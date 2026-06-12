"""
Multi-panel pathway convergence figure: Hallmark ORA on htNSC target union vs
mechanotransduction priors (Piezo1 / LIFU-motivated GMTs), heatmap of gene-set Jaccard overlap.

Supplementary UpSet: genes shared across multiple stress Hallmark terms (ORA gene lists).

Outputs:
  <outputs-dir>/figures/sa_bundle/fig_pathway_convergence.png
  <outputs-dir>/figures/sa_bundle/fig_supp_hallmark_stress_upset.png
  <outputs-dir>/exploratory_pathway_convergence_jaccard.csv
  <outputs-dir>/exploratory_pathway_convergence_hallmark_upset_intersections.csv
  <outputs-dir>/exploratory_pathway_convergence_hallmark_upset_membership.csv

Run from feasibility_study/:
  python figures/build_pathway_convergence_figure.py
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sa_nsc_lifu_computational_suite import CURATED_GENE_SETS
from journal_tier_string_bridge import MECHANISM_SEEDS

# Hallmark terms to highlight (must match MSigDB Hallmark names in ORA output).
STRESS_HALLMARK_TERMS: tuple[str, ...] = (
    "TGF-beta Signaling",
    "Hypoxia",
    "TNF-alpha Signaling via NF-kB",
    "Apoptosis",
    "Inflammatory Response",
)

HALLMARK_SHORT_LABELS: dict[str, str] = {
    "TGF-beta Signaling": "TGF-β",
    "Hypoxia": "Hypoxia",
    "TNF-alpha Signaling via NF-kB": "TNF/NF-κB",
    "Apoptosis": "Apoptosis",
    "Inflammatory Response": "Inflammatory",
}


def _norm(s: str) -> str:
    return str(s).strip().upper()


def _parse_gene_list(cell: object) -> set[str]:
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return set()
    return {_norm(g) for g in str(cell).split(";") if str(g).strip()}


def _load_hallmark_sets(ora: pd.DataFrame) -> dict[str, set[str]]:
    hallmark_sets: dict[str, set[str]] = {}
    for _, row in ora.iterrows():
        term = str(row.get("Term", "")).strip()
        if term not in STRESS_HALLMARK_TERMS:
            continue
        hallmark_sets[term] = _parse_gene_list(row.get("Genes"))
    return hallmark_sets


def _upset_intersections(
    hallmark_sets: dict[str, set[str]],
    *,
    min_terms: int = 2,
) -> tuple[list[str], list[tuple[tuple[bool, ...], list[str]]]]:
    """Return set names and (membership pattern, genes) for genes in >= min_terms Hallmark lists."""
    names = list(hallmark_sets.keys())
    buckets: dict[tuple[bool, ...], list[str]] = defaultdict(list)
    for g in sorted(set().union(*hallmark_sets.values())):
        pat = tuple(g in hallmark_sets[n] for n in names)
        if sum(pat) >= min_terms:
            buckets[pat].append(g)
    ranked = sorted(buckets.items(), key=lambda kv: len(kv[1]), reverse=True)
    return names, ranked


def _write_upset_tables(
    out_dir: Path,
    names: list[str],
    ranked: list[tuple[tuple[bool, ...], list[str]]],
) -> None:
    short = [HALLMARK_SHORT_LABELS.get(n, n) for n in names]
    rows_int: list[dict] = []
    rows_mem: list[dict] = []
    for idx, (pat, genes) in enumerate(ranked):
        active = [short[i] for i, a in enumerate(pat) if a]
        rows_int.append(
            {
                "intersection_id": idx,
                "n_hallmark_terms": len(active),
                "hallmark_terms": ";".join(active),
                "n_genes": len(genes),
                "genes": ";".join(genes[:500]),
            }
        )
        for g in genes:
            rows_mem.append(
                {
                    "gene": g,
                    "n_hallmark_terms": len(active),
                    "hallmark_terms": ";".join(active),
                    "intersection_id": idx,
                }
            )
    pd.DataFrame(rows_int).to_csv(
        out_dir / "exploratory_pathway_convergence_hallmark_upset_intersections.csv",
        index=False,
    )
    pd.DataFrame(rows_mem).to_csv(
        out_dir / "exploratory_pathway_convergence_hallmark_upset_membership.csv",
        index=False,
    )


def build_hallmark_stress_upset_figure(
    out_dir: Path,
    hallmark_sets: dict[str, set[str]],
    *,
    top_n: int = 14,
    min_terms: int = 2,
) -> bool:
    """Supplementary UpSet-style panel: multi-Hallmark gene intersections within ORA gene lists."""
    if len(hallmark_sets) < 2:
        return False
    names, ranked = _upset_intersections(hallmark_sets, min_terms=min_terms)
    if not ranked:
        return False
    _write_upset_tables(out_dir, names, ranked)
    ranked = ranked[:top_n]
    short = [HALLMARK_SHORT_LABELS.get(n, n) for n in names]
    n_sets = len(names)

    fig = plt.figure(figsize=(10.5, 5.6), dpi=150)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.08)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_bar)

    counts = [len(genes) for _, genes in ranked]
    x = np.arange(len(ranked))
    ax_bar.bar(x, counts, color="#1f4e79", alpha=0.9, width=0.72)
    ax_bar.set_ylabel("Genes in intersection")
    ax_bar.set_title(
        "UpSet: genes shared across ≥2 stress Hallmark terms (htNSC miRNA target union ORA gene lists)",
        fontsize=10,
        fontweight="bold",
    )
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    for j, (pat, _) in enumerate(ranked):
        active = [i for i, a in enumerate(pat) if a]
        for i in active:
            ax_mat.scatter(j, n_sets - 1 - i, s=58, c="#1f4e79", zorder=3)
        if len(active) >= 2:
            y0, y1 = n_sets - 1 - max(active), n_sets - 1 - min(active)
            ax_mat.plot([j, j], [y0, y1], color="#1f4e79", lw=2.2, zorder=2)

    ax_mat.set_yticks(np.arange(n_sets))
    ax_mat.set_yticklabels(list(reversed(short)), fontsize=9)
    ax_mat.set_xlabel("Intersection (sorted by size)")
    ax_mat.set_xlim(-0.6, len(ranked) - 0.4)
    ax_mat.spines["top"].set_visible(False)
    ax_mat.spines["right"].set_visible(False)
    plt.setp(ax_bar.get_xticklabels(), visible=False)

    out_png = out_dir / "figures" / "sa_bundle" / "fig_supp_hallmark_stress_upset.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")
    return True


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return float("nan")
    u = a | b
    if not u:
        return float("nan")
    return len(a & b) / len(u)


def build_pathway_convergence_figure(out_dir: Path) -> bool:
    out_dir = out_dir.resolve()
    ora_path = out_dir / "enrichr_hallmark_ora_mirtarbase_union.csv"
    union_path = out_dir / "mirna_target_union_genes.csv"
    string_path = out_dir / "exploratory_string_piezo1_bridge_summary.json"
    if not ora_path.is_file() or not union_path.is_file():
        return False

    ora = pd.read_csv(ora_path)
    union_genes: set[str] = set()
    udf = pd.read_csv(union_path)
    gcol = "gene" if "gene" in udf.columns else udf.columns[0]
    union_genes = {_norm(g) for g in udf[gcol].astype(str)}

    hallmark_sets = _load_hallmark_sets(ora)
    if not hallmark_sets:
        return False
    build_hallmark_stress_upset_figure(out_dir, hallmark_sets)

    mech_cols: dict[str, set[str]] = {
        "Piezo1_mechanism_seeds": {_norm(g) for g in MECHANISM_SEEDS},
        "Neuroinflammatory_signaling": {_norm(g) for g in CURATED_GENE_SETS.get("Neuroinflammatory_signaling", [])},
        "LIFU_mechanosensory_calcium": {_norm(g) for g in CURATED_GENE_SETS.get("LIFU_mechanosensory_calcium", [])},
        "miRNA_target_union": union_genes,
    }

    row_labels = list(hallmark_sets.keys())
    col_labels = list(mech_cols.keys())
    mat = np.zeros((len(row_labels), len(col_labels)), dtype=float)
    rows_out: list[dict] = []
    for i, rlab in enumerate(row_labels):
        rs = hallmark_sets[rlab]
        for j, clab in enumerate(col_labels):
            cs = mech_cols[clab]
            jac = _jaccard(rs, cs)
            mat[i, j] = jac if np.isfinite(jac) else 0.0
            rows_out.append(
                {
                    "hallmark_term": rlab,
                    "comparison_set": clab,
                    "jaccard": jac,
                    "n_hallmark_genes": len(rs),
                    "n_comparison_genes": len(cs),
                    "n_intersection": len(rs & cs),
                }
            )

    jaccard_df = pd.DataFrame(rows_out)
    jaccard_df.to_csv(out_dir / "exploratory_pathway_convergence_jaccard.csv", index=False)

    ora_sub = ora[ora["Term"].isin(STRESS_HALLMARK_TERMS)].copy()
    if ora_sub.empty:
        return False
    ora_sub["neglog10_padj"] = -np.log10(np.maximum(ora_sub["Adjusted P-value"].astype(float), 1e-300))
    ora_sub = ora_sub.sort_values("neglog10_padj", ascending=True)

    string_note = ""
    if string_path.is_file():
        try:
            sj = json.loads(string_path.read_text(encoding="utf-8"))
            n_ed = sj.get("n_edges_mechanism_seeds_to_union_subset")
            p_ed = sj.get("perm_p_edges_ge_obs_uniform_sample_null")
            if n_ed is not None:
                string_note = (
                    f"STRING (exploratory): {n_ed} edges (score≥filter) between "
                    f"Piezo1-mechanism seeds and high-burden union subset; "
                    f"uniform-size edge null P={p_ed}"
                )
        except Exception:
            pass

    fig = plt.figure(figsize=(12.5, 5.8), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.35], wspace=0.38)

    ax_a = fig.add_subplot(gs[0, 0])
    terms_a = ora_sub["Term"].astype(str)
    ypos = np.arange(len(terms_a))
    ax_a.barh(ypos, ora_sub["neglog10_padj"].values, color="#1f4e79", alpha=0.88, height=0.65)
    ax_a.set_yticks(ypos)
    ax_a.set_yticklabels(terms_a, fontsize=9)
    ax_a.set_xlabel(r"$-\log_{10}$(BH-adjusted P)")
    ax_a.set_title("A  Hallmark ORA on htNSC miRNA target union", fontsize=10, fontweight="bold")

    ax_b = fig.add_subplot(gs[0, 1])
    im = ax_b.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(0.15, float(np.nanmax(mat)) * 1.05))
    ax_b.set_xticks(np.arange(len(col_labels)))
    ax_b.set_xticklabels(
        [c.replace("_", " ") for c in col_labels],
        rotation=35,
        ha="right",
        fontsize=8,
    )
    ax_b.set_yticks(np.arange(len(row_labels)))
    ax_b.set_yticklabels(row_labels, fontsize=9)
    ax_b.set_title("B  Gene-set Jaccard (Hallmark genes vs priors / union)", fontsize=10, fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax_b.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8, color="black")
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard index")

    if string_note:
        fig.text(0.5, 0.02, string_note, ha="center", fontsize=8, style="italic", wrap=True)

    fig.suptitle(
        "Pathway convergence: htNSC targetome stress axes and Piezo1–LIFU mechanotransduction neighborhood",
        fontsize=11,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0.06 if string_note else 0, 1, 0.98])

    out_png = out_dir / "figures" / "sa_bundle" / "fig_pathway_convergence.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_dir / 'exploratory_pathway_convergence_jaccard.csv'}")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    args = ap.parse_args()
    ok = build_pathway_convergence_figure(args.outputs_dir)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
