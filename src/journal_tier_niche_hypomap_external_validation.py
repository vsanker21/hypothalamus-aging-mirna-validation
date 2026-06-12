"""
External third-ventricle niche validation via HypoMap (Steuernagel et al., Nat Metab 2022).

Independent of internal marker-module labels: uses Spearman correlations between
GSE188646 cluster mean expression and HypoMap C185 named cell-type means
(outputs/gse188646_hypomap_mapping/hypomap_custom_ref_spearman.csv from
r/gse188646_hypomap_reference_mapping.R + data/references/hypomap_cellxgene_C185_named_mean_X_min200.csv).

Third-ventricle / circumventricular reference types (curated from HypoMap taxonomy):
  *Tanycytes, *Ependymal, *ParsTuber (median eminence / third-ventricle floor).

Tests:
  1) Per-stratum external niche rho and best-matching HypoMap type
  2) Concordance with internal marker-module niche panel (Jaccard, Fisher enrichment)
  3) |logFC| attenuation (delta_median_abs_logfc) in HypoMap-validated vs other strata
  4) Concordant-only stratum analysis (marker niche AND external rho >= threshold)

Outputs:
  exploratory_niche_hypomap_external_validation_per_stratum.csv
  exploratory_niche_hypomap_external_validation_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# HypoMap C185 types annotating third-ventricle / ME niche (Steuernagel HypoMap harmonized taxonomy).
HYPOMAP_NICHE_TYPE_RE = r"Tanycyte|Ependymal|ParsTuber"

DEFAULT_RHO_THRESHOLD = 0.82


def _fisher_enrichment(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    try:
        tab = [[a, b], [c, d]]
        or_, p = stats.fisher_exact(tab, alternative="two-sided")
        return float(or_), float(p)
    except Exception:
        return float("nan"), float("nan")


def _wilcoxon_two_sample(a: pd.Series, b: pd.Series) -> float:
    a = a.dropna().astype(float)
    b = b.dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def run_niche_hypomap_external_validation(
    out_dir: Path,
    log,
    *,
    rho_threshold: float = DEFAULT_RHO_THRESHOLD,
    n_perm: int = 2000,
    seed: int = 52,
) -> None:
    log("\n=== External HypoMap third-ventricle niche validation (independent reference) ===")
    hm_path = out_dir / "gse188646_hypomap_mapping" / "hypomap_custom_ref_spearman.csv"
    niche_path = out_dir / "exploratory_niche_lability_per_stratum.csv"
    mod_path = out_dir / "gse188646_cluster_annotation" / "cluster_module_scores.csv"

    if not hm_path.is_file():
        log(
            "HypoMap external validation: missing hypomap_custom_ref_spearman.csv — "
            "run r/gse188646_hypomap_reference_mapping.R with GSE188646_HYPOMAP_REF_EXPR_CSV; skipped."
        )
        return
    if not niche_path.is_file():
        log("HypoMap external validation: missing exploratory_niche_lability_per_stratum.csv; skipped.")
        return

    hm = pd.read_csv(hm_path)
    hm["seurat_cluster_id"] = hm["seurat_cluster_id"].astype(str).str.strip()
    niche_df = pd.read_csv(niche_path)
    niche_df["stratum"] = niche_df["stratum"].astype(str).str.strip()

    hm_niche = hm[hm["ref_celltype"].astype(str).str.contains(HYPOMAP_NICHE_TYPE_RE, case=False, na=False, regex=True)].copy()
    if hm_niche.empty:
        log("HypoMap external validation: no third-ventricle reference types matched; skipped.")
        return

    ext_rows: list[dict] = []
    for sid, g in hm.groupby("seurat_cluster_id"):
        g_all = g.sort_values("rho", ascending=False)
        g_n = hm_niche[hm_niche["seurat_cluster_id"] == sid].sort_values("rho", ascending=False)
        if g_n.empty:
            continue
        best_n = g_n.iloc[0]
        rank1_all = g_all.iloc[0]
        ext_rows.append(
            {
                "stratum": sid,
                "hypomap_niche_best_type": best_n["ref_celltype"],
                "hypomap_niche_best_rho": float(best_n["rho"]),
                "hypomap_niche_best_padj": float(best_n.get("p_adj", np.nan)),
                "hypomap_rank1_all_types": rank1_all["ref_celltype"],
                "hypomap_rank1_all_rho": float(rank1_all["rho"]),
                "hypomap_niche_in_top3_all_types": bool(
                    best_n["ref_celltype"] in g_all.head(3)["ref_celltype"].values
                ),
            }
        )

    ext = pd.DataFrame(ext_rows)
    d = niche_df.merge(ext, on="stratum", how="left")
    d["hypomap_niche_validated"] = d["hypomap_niche_best_rho"].astype(float) >= rho_threshold
    d["concordant_marker_and_hypomap"] = (
        d["is_third_ventricle_niche"].astype(bool) & d["hypomap_niche_validated"].astype(bool)
    )

    if mod_path.is_file():
        mod = pd.read_csv(mod_path)
        idcol = "seurat_cluster_id" if "seurat_cluster_id" in mod.columns else "cluster"
        mod[idcol] = mod[idcol].astype(str).str.strip()
        tv_cols = [c for c in ("Tanycyte_ependymal", "NSC_like", "Radial_glia_like") if c in mod.columns]
        if tv_cols:
            mod["third_ventricle_module_score"] = mod[tv_cols].astype(float).mean(axis=1)
            d = d.merge(
                mod[[idcol, "third_ventricle_module_score"]].rename(columns={idcol: "stratum"}),
                on="stratum",
                how="left",
            )

    out_csv = out_dir / "exploratory_niche_hypomap_external_validation_per_stratum.csv"
    d.to_csv(out_csv, index=False)

    n_de = int(d["delta_median_abs_logfc"].notna().sum())
    n_marker = int(d["is_third_ventricle_niche"].astype(bool).sum())
    n_hypomap = int(d["hypomap_niche_validated"].astype(bool).sum())
    n_concordant = int(d["concordant_marker_and_hypomap"].sum())
    marker_ids = set(d.loc[d["is_third_ventricle_niche"].astype(bool), "stratum"])
    hypomap_ids = set(d.loc[d["hypomap_niche_validated"].astype(bool), "stratum"])
    concordant_ids = sorted(marker_ids & hypomap_ids, key=lambda x: int(x) if x.isdigit() else x)
    jaccard = (
        len(marker_ids & hypomap_ids) / len(marker_ids | hypomap_ids)
        if (marker_ids or hypomap_ids)
        else float("nan")
    )

    # Fisher: marker-niche vs hypomap-validated among DE strata
    a = int(((d["is_third_ventricle_niche"].astype(bool)) & (d["hypomap_niche_validated"].astype(bool))).sum())
    b = int((d["is_third_ventricle_niche"].astype(bool) & ~d["hypomap_niche_validated"].astype(bool)).sum())
    c = int((~d["is_third_ventricle_niche"].astype(bool) & d["hypomap_niche_validated"].astype(bool)).sum())
    d_ct = int((~d["is_third_ventricle_niche"].astype(bool) & ~d["hypomap_niche_validated"].astype(bool)).sum())
    fisher_or, fisher_p = _fisher_enrichment(a, b, c, d_ct)
    if fisher_or is not None and not np.isfinite(fisher_or):
        fisher_or = None

    hv = d[d["hypomap_niche_validated"].astype(bool)]
    other = d[~d["hypomap_niche_validated"].astype(bool)]
    wilcox_hypomap = _wilcoxon_two_sample(
        hv["delta_median_abs_logfc"], other["delta_median_abs_logfc"]
    )

    conc = d[d["concordant_marker_and_hypomap"]]
    non_conc = d[~d["concordant_marker_and_hypomap"].astype(bool)]
    wilcox_concordant = _wilcoxon_two_sample(
        conc["delta_median_abs_logfc"], non_conc["delta_median_abs_logfc"]
    )

    # Permutation: shuffle hypomap_niche_validated labels
    rng = np.random.default_rng(seed)
    labels = d["hypomap_niche_validated"].astype(bool).values
    deltas = d["delta_median_abs_logfc"].astype(float).values
    obs = (
        float(np.nanmedian(deltas[labels]) - np.nanmedian(deltas[~labels]))
        if labels.any() and (~labels).any()
        else float("nan")
    )
    perm_diffs = []
    for _ in range(n_perm):
        lp = rng.permutation(labels)
        if lp.sum() < 1 or (~lp).sum() < 1:
            continue
        perm_diffs.append(float(np.nanmedian(deltas[lp]) - np.nanmedian(deltas[~lp])))
    perm_diffs = np.asarray(perm_diffs, dtype=float)
    perm_p = float(np.mean(perm_diffs >= obs)) if len(perm_diffs) and np.isfinite(obs) else float("nan")

    # Spearman: external niche rho vs attenuation relief (higher delta = less attenuation)
    rho_sp, p_sp = stats.spearmanr(
        d["hypomap_niche_best_rho"].astype(float),
        d["delta_median_abs_logfc"].astype(float),
        nan_policy="omit",
    )

    summ = {
        "reference_atlas": "HypoMap CELLxGENE C185_named (Steuernagel et al., Nat Metab 2022)",
        "reference_matrix": "data/references/hypomap_cellxgene_C185_named_mean_X_min200.csv",
        "hypomap_mapping_source": "outputs/gse188646_hypomap_mapping/hypomap_custom_ref_spearman.csv",
        "niche_type_pattern": "Tanycytes|Ependymal|ParsTuber",
        "rho_threshold_hypomap_validated": rho_threshold,
        "n_strata_with_de": n_de,
        "n_marker_niche_strata": n_marker,
        "n_hypomap_validated_strata": n_hypomap,
        "n_concordant_strata": n_concordant,
        "concordant_strata_ids": concordant_ids,
        "marker_niche_strata_ids": sorted(marker_ids, key=lambda x: int(x) if str(x).isdigit() else str(x)),
        "hypomap_validated_strata_ids": sorted(hypomap_ids, key=lambda x: int(x) if str(x).isdigit() else str(x)),
        "jaccard_marker_vs_hypomap_niche_ids": float(jaccard),
        "fisher_marker_and_hypomap_a": a,
        "fisher_marker_only_b": b,
        "fisher_hypomap_only_c": c,
        "fisher_neither_d": d_ct,
        "fisher_odds_ratio_marker_vs_hypomap": fisher_or,
        "fisher_p_marker_vs_hypomap": fisher_p,
        "cluster_22_hypomap_best_type": d.loc[d["stratum"] == "22", "hypomap_niche_best_type"].iloc[0]
        if (d["stratum"] == "22").any()
        else None,
        "cluster_22_hypomap_niche_rho": float(d.loc[d["stratum"] == "22", "hypomap_niche_best_rho"].iloc[0])
        if (d["stratum"] == "22").any()
        else None,
        "cluster_6_hypomap_niche_best_type": d.loc[d["stratum"] == "6", "hypomap_niche_best_type"].iloc[0]
        if (d["stratum"] == "6").any()
        else None,
        "cluster_6_hypomap_niche_rho": float(d.loc[d["stratum"] == "6", "hypomap_niche_best_rho"].iloc[0])
        if (d["stratum"] == "6").any()
        else None,
        "wilcoxon_p_delta_median_hypomap_validated_vs_other": wilcox_hypomap,
        "wilcoxon_p_delta_median_concordant_vs_rest": wilcox_concordant,
        "perm_n_hypomap_label_shuffle": n_perm,
        "perm_p_median_delta_diff_hypomap_ge_obs": perm_p,
        "spearman_rho_hypomap_niche_score_vs_delta_median_abs_logfc": float(rho_sp),
        "spearman_p_hypomap_niche_score_vs_delta_median_abs_logfc": float(p_sp),
        "methodology_note": (
            "External validation uses an independent hypothalamus atlas (HypoMap; different animals, "
            "protocols, and normalization than GSE188646). Spearman rho compares GSE188646 cluster "
            "mean expression to HypoMap cell-type mean profiles — not label transfer. "
            "Third-ventricle reference types are HypoMap Tanycytes, Ependymal, and ParsTuber clusters. "
            "Concordant strata satisfy BOTH internal marker-module niche rules AND rho >= threshold to "
            "the external reference."
        ),
        "caveat": (
            "HypoMap is not spatially resolved in situ; ParsTuber types are median-eminence-adjacent, "
            "not a direct third-ventricle dissection. High global correlation between many hypothalamic "
            "clusters and Tanycyte profiles can occur because of shared glial/ventricular gene programs — "
            "use concordance with marker modules and per-cluster biology, not rho alone."
        ),
    }
    json_path = out_dir / "exploratory_niche_hypomap_external_validation_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {out_csv.name} ({len(d)} strata).")
    log(f"Wrote {json_path.name}")
    or_txt = f"{fisher_or:.3g}" if fisher_or is not None and np.isfinite(fisher_or) else "inf"
    log(
        f"Marker niche n={n_marker}; HypoMap-validated (rho>={rho_threshold}) n={n_hypomap}; "
        f"concordant n={n_concordant} (ids={concordant_ids}); Fisher OR={or_txt}, P={fisher_p:.3g}"
    )
    if summ.get("cluster_22_hypomap_best_type"):
        log(
            f"Cluster 22 external match: {summ['cluster_22_hypomap_best_type']} "
            f"(rho={summ['cluster_22_hypomap_niche_rho']:.3f})"
        )
