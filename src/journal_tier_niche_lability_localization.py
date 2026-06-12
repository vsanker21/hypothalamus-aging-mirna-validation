"""
Third-ventricle / tanycyte–NSC niche localization of miRNA-target aging lability.

Uses per-cluster (not whole-tissue) pseudobulk DE from outputs/gse188646_strata/* joined to
marker-module labels (outputs/gse188646_cluster_annotation/cluster_putative_labels.csv).

For each stratum with DE:
  - Mann–Whitney |logFC| in miRNA-target union vs non-union (same test as global cross-modal)
  - Effect size: delta median |logFC| (union − non-union) and rank-biserial r

Cross-stratum inference (exploratory):
  - Compare lability effect sizes between third-ventricle niche strata vs all other strata
    (Wilcoxon on delta_median_abs_logfc; cell-count–weighted variant)
  - Permutation null: shuffle niche labels across strata

Outputs:
  exploratory_niche_lability_per_stratum.csv
  exploratory_niche_lability_localization_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

NICHE_MODULES = frozenset(
    {
        "Tanycyte_ependymal",
        "NSC_like",
        "Radial_glia_like",
    }
)
ASTRO_NICHE_PAIR = frozenset({"Astrocyte"})


def _assign_niche_class(row: pd.Series) -> str:
    r1 = str(row.get("rank1_module") or "").strip()
    r2 = str(row.get("rank2_module") or "").strip()
    r1z = float(row.get("rank1_z")) if pd.notna(row.get("rank1_z")) else float("nan")
    r2z = float(row.get("rank2_z")) if pd.notna(row.get("rank2_z")) else float("nan")

    if r1 in NICHE_MODULES:
        return "third_ventricle_niche"
    if r2 in NICHE_MODULES and np.isfinite(r2z) and r2z >= 0.8:
        if r1 in ASTRO_NICHE_PAIR or (np.isfinite(r1z) and r2z >= r1z - 1.2):
            return "third_ventricle_niche"
    if r1 in ASTRO_NICHE_PAIR and r2 in NICHE_MODULES:
        return "astrocyte_adjacent_niche"
    return "other_hypothalamic"


def _stratum_lability_stats(deg: pd.DataFrame, union_set: set[str]) -> dict:
    deg = deg.copy()
    deg["gene_u"] = deg["gene"].astype(str).str.strip().str.upper()
    in_u = deg["gene_u"].isin(union_set)
    u_fc = np.abs(deg.loc[in_u, "logFC"].astype(float).values)
    o_fc = np.abs(deg.loc[~in_u, "logFC"].astype(float).values)
    if len(u_fc) < 5 or len(o_fc) < 5:
        return {}
    med_u, med_o = float(np.median(u_fc)), float(np.median(o_fc))
    try:
        mw = stats.mannwhitneyu(u_fc, o_fc, alternative="two-sided")
        mw_p = float(mw.pvalue)
        # Rank-biserial correlation (direction: union more labile if positive)
        n1, n2 = len(u_fc), len(o_fc)
        r_rb = 1.0 - (2.0 * mw.statistic) / (n1 * n2)
    except ValueError:
        mw_p, r_rb = float("nan"), float("nan")
    return {
        "n_genes_union": int(in_u.sum()),
        "n_genes_nonunion": int((~in_u).sum()),
        "median_abs_logfc_union": med_u,
        "median_abs_logfc_nonunion": med_o,
        "delta_median_abs_logfc": med_u - med_o,
        "mannwhitney_abs_logFC_union_vs_nonunion_p": mw_p,
        "rank_biserial_abs_logfc": r_rb,
    }


def run_niche_lability_localization(
    out_dir: Path,
    log,
    *,
    n_perm: int = 2000,
    seed: int = 51,
) -> None:
    log("\n=== Third-ventricle niche localization of miRNA-target aging lability (per-cluster DE) ===")
    summ_path = out_dir / "exploratory_crossmodal_celltype_strata_summary.csv"
    labels_path = out_dir / "gse188646_cluster_annotation" / "cluster_putative_labels.csv"
    long_path = out_dir / "mirna_targets_long.csv"
    strata_dir = out_dir / "gse188646_strata"

    if not summ_path.is_file():
        log("Niche lability: missing exploratory_crossmodal_celltype_strata_summary.csv; skipped.")
        return
    if not labels_path.is_file():
        log("Niche lability: missing cluster_putative_labels.csv — run cluster marker mapping R script; skipped.")
        return
    if not long_path.is_file():
        log("Niche lability: missing mirna_targets_long.csv; skipped.")
        return

    long_df = pd.read_csv(long_path)
    tcol = "target_gene" if "target_gene" in long_df.columns else long_df.columns[-1]
    union_set = {str(g).strip().upper() for g in long_df[tcol].astype(str)}

    summ = pd.read_csv(summ_path)
    summ["stratum"] = summ["stratum"].astype(str).str.strip()
    lab = pd.read_csv(labels_path)
    idcol = "seurat_cluster_id" if "seurat_cluster_id" in lab.columns else "stratum"
    lab[idcol] = lab[idcol].astype(str).str.strip()
    label_cols = [
        c
        for c in (
            "rank1_module",
            "rank2_module",
            "rank1_z",
            "rank2_z",
            "delta_z",
            "cluster",
        )
        if c in lab.columns
    ]
    d = summ.merge(
        lab[[idcol] + label_cols].rename(columns={idcol: "stratum"}),
        on="stratum",
        how="left",
        suffixes=("", "_lab"),
    )

    rows: list[dict] = []
    for _, row in d.iterrows():
        fn = str(row.get("deg_csv") or "").strip()
        if not fn:
            continue
        deg_path = strata_dir / fn
        if not deg_path.is_file():
            continue
        deg = pd.read_csv(deg_path)
        if "logFC" not in deg.columns:
            continue
        st = _stratum_lability_stats(deg, union_set)
        if not st:
            continue
        niche = _assign_niche_class(row)
        rows.append(
            {
                "stratum": row["stratum"],
                "stratum_col": row.get("stratum_col", ""),
                "rank1_module": row.get("rank1_module"),
                "rank2_module": row.get("rank2_module"),
                "rank1_z": row.get("rank1_z"),
                "rank2_z": row.get("rank2_z"),
                "delta_z": row.get("delta_z"),
                "n_cells": row.get("n_cells"),
                "niche_class": niche,
                "is_third_ventricle_niche": niche in ("third_ventricle_niche", "astrocyte_adjacent_niche"),
                "deg_csv": fn,
                **st,
            }
        )

    if not rows:
        log("Niche lability: no strata with computable lability stats; skipped.")
        return

    df = pd.DataFrame(rows)
    out_csv = out_dir / "exploratory_niche_lability_per_stratum.csv"
    df.to_csv(out_csv, index=False)

    niche = df[df["is_third_ventricle_niche"]]
    other = df[~df["is_third_ventricle_niche"]]
    global_mw_p = None
    cm_path = out_dir / "exploratory_crossmodal_mirna_aging_summary.json"
    if cm_path.is_file():
        try:
            global_mw_p = json.loads(cm_path.read_text(encoding="utf-8")).get(
                "mannwhitney_abs_logFC_union_vs_nonunion_p"
            )
        except Exception:
            pass

    def _wilcox(a: pd.Series, b: pd.Series) -> tuple[float, float]:
        a = a.dropna().astype(float)
        b = b.dropna().astype(float)
        if len(a) < 2 or len(b) < 2:
            return float("nan"), float("nan")
        w = stats.mannwhitneyu(a, b, alternative="two-sided")
        return float(w.statistic), float(w.pvalue)

    u_stat, wilcox_p_delta = _wilcox(niche["delta_median_abs_logfc"], other["delta_median_abs_logfc"])
    _, wilcox_p_rb = _wilcox(niche["rank_biserial_abs_logfc"], other["rank_biserial_abs_logfc"])

    # Permutation: shuffle is_third_ventricle_niche labels
    rng = np.random.default_rng(seed)
    labels = df["is_third_ventricle_niche"].astype(bool).values
    deltas = df["delta_median_abs_logfc"].astype(float).values
    obs_diff = float(np.nanmedian(deltas[labels]) - np.nanmedian(deltas[~labels])) if labels.any() and (~labels).any() else float("nan")
    perm_diffs = []
    for _ in range(n_perm):
        lp = rng.permutation(labels)
        if lp.sum() < 1 or (~lp).sum() < 1:
            continue
        perm_diffs.append(float(np.nanmedian(deltas[lp]) - np.nanmedian(deltas[~lp])))
    perm_diffs = np.asarray(perm_diffs, dtype=float)
    if len(perm_diffs) and np.isfinite(obs_diff):
        perm_p = float(np.mean(perm_diffs >= obs_diff))
    else:
        perm_p = float("nan")

    # Cell-weighted mean delta (niche strata contribute proportionally to nucleus count)
    w_n = df.loc[df["is_third_ventricle_niche"], "n_cells"].astype(float)
    w_o = df.loc[~df["is_third_ventricle_niche"], "n_cells"].astype(float)
    d_n = df.loc[df["is_third_ventricle_niche"], "delta_median_abs_logfc"].astype(float)
    d_o = df.loc[~df["is_third_ventricle_niche"], "delta_median_abs_logfc"].astype(float)
    wmean_n = float(np.average(d_n, weights=w_n)) if len(d_n) and w_n.sum() > 0 else float("nan")
    wmean_o = float(np.average(d_o, weights=w_o)) if len(d_o) and w_o.sum() > 0 else float("nan")

    summ_out = {
        "n_strata_total": int(len(df)),
        "n_third_ventricle_niche_strata": int(niche.shape[0]),
        "n_other_strata": int(other.shape[0]),
        "third_ventricle_niche_strata_ids": niche["stratum"].astype(str).tolist(),
        "median_delta_abs_logfc_niche": float(np.nanmedian(niche["delta_median_abs_logfc"])) if len(niche) else None,
        "median_delta_abs_logfc_other": float(np.nanmedian(other["delta_median_abs_logfc"])) if len(other) else None,
        "weighted_mean_delta_abs_logfc_niche": wmean_n,
        "weighted_mean_delta_abs_logfc_other": wmean_o,
        "obs_median_delta_diff_niche_minus_other": obs_diff,
        "wilcoxon_p_delta_median_niche_vs_other": wilcox_p_delta,
        "wilcoxon_p_rank_biserial_niche_vs_other": wilcox_p_rb,
        "perm_n_niche_label_shuffle": n_perm,
        "perm_p_median_delta_diff_ge_obs": perm_p,
        "global_pseudobulk_mannwhitney_p_reference": global_mw_p,
        "methodology_note": (
            "Per-cluster pseudobulk young vs aged DE (GSE188646 strata); lability = Mann–Whitney on "
            "|logFC| for miRTarBase union vs non-union within that cluster's gene universe. "
            "Third-ventricle niche strata: rank1 Tanycyte_ependymal/NSC_like/Radial_glia_like, or "
            "Astrocyte rank1 with Tanycyte/NSC/Radial rank2 (module z ≥ 0.8), or Astrocyte+Tanycyte pair. "
            "Cross-stratum Wilcoxon compares delta_median_abs_logfc between niche-labelled vs other clusters."
        ),
        "caveat": (
            "Clusters are Seurat partitions with heuristic marker labels, not spatial third-ventricle dissection; "
            "niche strata with few DE genes may be unstable. Exploratory only."
        ),
    }
    json_path = out_dir / "exploratory_niche_lability_localization_summary.json"
    json_path.write_text(json.dumps(summ_out, indent=2), encoding="utf-8")
    log(f"Wrote {out_csv.name} ({len(df)} strata).")
    log(f"Wrote {json_path.name}")
    log(
        f"Third-ventricle niche strata n={summ_out['n_third_ventricle_niche_strata']} "
        f"(ids={summ_out['third_ventricle_niche_strata_ids']}); "
        f"median d|logFC| niche={summ_out['median_delta_abs_logfc_niche']}, other={summ_out['median_delta_abs_logfc_other']}; "
        f"Wilcoxon P={wilcox_p_delta}; perm P={perm_p}"
    )
