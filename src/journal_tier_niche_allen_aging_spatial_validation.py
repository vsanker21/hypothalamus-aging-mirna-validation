"""
Jin et al. 2025 Resolve RSTE3 spatial transcriptomics (Brain Image Library) — V3 aging validation.

Downloads processed RSTE3 tables from BIL doi:10.35077/g.1157 and computes cell-level
age x V3 niche statistics (tanycytes / ependymal / astroependymal at hypothalamus ROI).

Outputs:
  exploratory_allen_rstE3_access_status.json
  exploratory_allen_aging_spatial_rstE3_gene_overlap.csv
  exploratory_allen_aging_spatial_rstE3_markers_in_s22.csv
  exploratory_allen_rstE3_age_by_subclass_summary.csv
  exploratory_allen_rstE3_cell_level_age_stats.csv
  exploratory_allen_aging_spatial_validation_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from bil_rstE3_fetch import BIL_DOI, BIL_DOWNLOAD_BASE, ensure_rstE3_files, read_rstE3_metadata

# RSTE3 panel genes highlighted in Jin et al. Fig. 4 / Ext Data Fig. 11 (V3 tanycytes/ependymal)
RSTE3_CURATED_MARKERS = {
    "OASL2", "IFIT1", "CCND2", "CTNNA2", "GPR50", "TM4SF1", "H2-K1", "IFI27",
    "FOXJ1", "RAX", "TTR", "FABP7", "DIO2", "PROM1", "VIM", "CLU",
}

V3_NICHE_SUBCLASSES = ("Tanycyte NN", "Ependymal NN", "Astroependymal NN")
AGE_YOUNG_DAYS = 56
AGE_OLD_DAYS = 540


def _try_bil_access(log) -> dict:
    status = {
        "doi": BIL_DOI,
        "download_base": BIL_DOWNLOAD_BASE,
        "accessible": False,
        "processed_files": list(ensure_rstE3_files(log=log).keys()),
        "notes": [],
    }
    try:
        paths = ensure_rstE3_files(log=log)
        meta = paths["metadata"]
        expr = paths["cellxgene"]
        if meta.is_file() and meta.stat().st_size > 1_000_000 and expr.is_file() and expr.stat().st_size > 1_000_000:
            status["accessible"] = True
            status["metadata_bytes"] = meta.stat().st_size
            status["cellxgene_bytes"] = expr.stat().st_size
            status["notes"].append("RSTE3 processed metadata + cellxgene cached locally.")
        else:
            status["notes"].append("BIL RSTE3 files not fully cached.")
    except Exception as exc:
        status["notes"].append(f"RSTE3 fetch failed: {exc}")
    return status


def _load_niche_de_genes(out_dir: Path, strata: list[str]) -> set[str]:
    genes: set[str] = set()
    padj_cols = ("FDR", "padj", "p_val_adj", "adj.P.Val")
    for sid in strata:
        p = out_dir / "gse188646_strata" / f"stratum_{sid}_young_vs_aged_deg.csv"
        if not p.is_file():
            continue
        deg = pd.read_csv(p)
        padj_col = next((c for c in padj_cols if c in deg.columns), None)
        for _, row in deg.iterrows():
            g = str(row.get("gene", "")).strip().upper()
            if not g:
                continue
            if padj_col and float(row[padj_col]) <= 0.05:
                genes.add(g)
    return genes


def _load_s22_marker_logfc(out_dir: Path) -> dict[str, float]:
    s22_path = out_dir / "gse188646_strata" / "stratum_22_young_vs_aged_deg.csv"
    if not s22_path.is_file():
        return {}
    deg = pd.read_csv(s22_path)
    deg["gene_u"] = deg["gene"].astype(str).str.upper()
    out: dict[str, float] = {}
    for _, row in deg.iterrows():
        lfc = row.get("logFC", row.get("log2FoldChange", np.nan))
        if np.isfinite(lfc):
            out[str(row["gene_u"])] = float(lfc)
    return out


def _age_label(days: int | float) -> str:
    d = int(days)
    if d == AGE_YOUNG_DAYS:
        return "young_2mo"
    if d == AGE_OLD_DAYS:
        return "aged_18mo"
    return f"age_{d}d"


def _compute_age_by_subclass_summary(meta_v3: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subclass, g in meta_v3.groupby("subclass_label"):
        for age_days, ga in g.groupby("age"):
            rows.append(
                {
                    "subclass_label": subclass,
                    "age_days": int(age_days),
                    "age_group": _age_label(age_days),
                    "n_cells": len(ga),
                    "mean_density": float(ga["density"].mean()) if "density" in ga.columns else np.nan,
                    "median_n_transcripts": float(ga["n_transcripts"].median())
                    if "n_transcripts" in ga.columns
                    else np.nan,
                    "mean_x": float(ga["x"].mean()) if "x" in ga.columns else np.nan,
                    "mean_y": float(ga["y"].mean()) if "y" in ga.columns else np.nan,
                }
            )
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary

    # Proportions within subclass (young vs aged)
    prop_rows = []
    for subclass in summary["subclass_label"].unique():
        sub = summary[summary["subclass_label"] == subclass]
        total = sub["n_cells"].sum()
        young_n = sub.loc[sub["age_days"] == AGE_YOUNG_DAYS, "n_cells"].sum()
        old_n = sub.loc[sub["age_days"] == AGE_OLD_DAYS, "n_cells"].sum()
        prop_rows.append(
            {
                "subclass_label": subclass,
                "n_cells_total": int(total),
                "n_young": int(young_n),
                "n_aged": int(old_n),
                "frac_young": float(young_n / total) if total else np.nan,
                "frac_aged": float(old_n / total) if total else np.nan,
            }
        )
    prop_df = pd.DataFrame(prop_rows)
    return summary.merge(prop_df, on="subclass_label", how="left")


def _mannwhitney_gene_tests(
    meta_v3: pd.DataFrame,
    expr_markers: pd.DataFrame,
    s22_logfc: dict[str, float],
) -> pd.DataFrame:
    """Per-gene Mann-Whitney (aged vs young) within each V3 subclass."""
    expr_t = expr_markers.T
    expr_t.index = expr_t.index.astype(str)
    merged = meta_v3.set_index("sample_id").join(expr_t, how="inner")
    gene_cols = [c for c in expr_t.columns if c in merged.columns]

    rows = []
    for subclass, g in merged.groupby("subclass_label"):
        young = g[g["age"] == AGE_YOUNG_DAYS]
        old = g[g["age"] == AGE_OLD_DAYS]
        for gene in gene_cols:
            yv = young[gene].astype(float)
            ov = old[gene].astype(float)
            if len(yv) < 5 or len(ov) < 5:
                continue
            try:
                mw = stats.mannwhitneyu(ov, yv, alternative="two-sided")
                p = float(mw.pvalue)
            except Exception:
                p = float("nan")
            mean_y = float(yv.mean())
            mean_o = float(ov.mean())
            pseudo = 0.5
            log2fc = float(np.log2((mean_o + pseudo) / (mean_y + pseudo)))
            gene_u = str(gene).upper()
            s22 = s22_logfc.get(gene_u, np.nan)
            direction_match = (
                np.sign(log2fc) == np.sign(s22)
                if np.isfinite(log2fc) and np.isfinite(s22) and log2fc != 0 and s22 != 0
                else np.nan
            )
            rows.append(
                {
                    "subclass_label": subclass,
                    "gene": gene,
                    "n_young": len(yv),
                    "n_aged": len(ov),
                    "mean_expr_young": mean_y,
                    "mean_expr_aged": mean_o,
                    "log2fc_aged_vs_young": log2fc,
                    "mannwhitney_p": p,
                    "gse188646_s22_logFC": s22,
                    "direction_concordant_with_s22": direction_match,
                    "pct_expressing_young": float((yv > 0).mean()),
                    "pct_expressing_aged": float((ov > 0).mean()),
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    mask = out["mannwhitney_p"].notna()
    if mask.any():
        _, q, _, _ = multipletests(out.loc[mask, "mannwhitney_p"].astype(float), method="fdr_bh")
        out.loc[mask, "fdr_bh"] = q
    return out.sort_values(["subclass_label", "mannwhitney_p"])


def _records_json_safe(df: pd.DataFrame) -> list[dict]:
    out: list[dict] = []
    for rec in df.to_dict(orient="records"):
        clean = {}
        for k, v in rec.items():
            if isinstance(v, (float, np.floating)) and not np.isfinite(v):
                clean[k] = None
            elif isinstance(v, (np.bool_, bool)):
                clean[k] = bool(v)
            elif isinstance(v, (np.integer, int)):
                clean[k] = int(v)
            else:
                clean[k] = v
        out.append(clean)
    return out


def _v3_spatial_interface_mask(meta_v3: pd.DataFrame) -> pd.Series:
    """
    Ventricular-interface mask: per RSTE3 slide (section), y between ependymal Q75 and
    tanycyte Q25; fallback to global HY niche quantiles when a slide lacks both types.
    """
    meta = meta_v3.copy()
    meta["_slide_id"] = meta["sample_id"].astype(str).map(_rstE3_slide_id)
    mask = pd.Series(False, index=meta_v3.index)

    epen_all = meta.loc[meta["subclass_label"] == "Ependymal NN", "y"].astype(float)
    tany_all = meta.loc[meta["subclass_label"] == "Tanycyte NN", "y"].astype(float)
    global_lo = global_hi = None
    if len(epen_all) >= 5 and len(tany_all) >= 5:
        global_lo = float(min(epen_all.quantile(0.75), tany_all.quantile(0.25)))
        global_hi = float(max(epen_all.quantile(0.75), tany_all.quantile(0.25)))
        if global_lo > global_hi:
            global_lo, global_hi = global_hi, global_lo

    for slide_id, g in meta.groupby("_slide_id"):
        epen = g[g["subclass_label"] == "Ependymal NN"]["y"].astype(float)
        tany = g[g["subclass_label"] == "Tanycyte NN"]["y"].astype(float)
        if len(epen) >= 5 and len(tany) >= 5:
            y_lo = float(min(epen.quantile(0.75), tany.quantile(0.25)))
            y_hi = float(max(epen.quantile(0.75), tany.quantile(0.25)))
        elif global_lo is not None:
            y_lo, y_hi = global_lo, global_hi
        else:
            continue
        if y_lo > y_hi:
            y_lo, y_hi = y_hi, y_lo
        pad = 0.05 * (y_hi - y_lo + 1.0)
        in_band = (g["y"].astype(float) >= y_lo - pad) & (g["y"].astype(float) <= y_hi + pad)
        mask.loc[g.index] = in_band
    return mask


def _rstE3_slide_id(sample_id: str) -> str:
    parts = str(sample_id).split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return str(sample_id)


def _compute_rstE3_cell_level_stats(out_dir: Path, log) -> dict | None:
    try:
        meta = read_rstE3_metadata(log=log)
    except Exception as exc:
        log(f"  RSTE3 metadata load failed: {exc}")
        return None

    v3_mask = (
        meta["roi"].astype(str).str.lower().eq("hy")
        & meta["subclass_label"].astype(str).isin(V3_NICHE_SUBCLASSES)
    )
    meta_v3 = meta.loc[v3_mask].copy()
    if meta_v3.empty:
        log("  No V3 niche cells after ROI/subclass filter.")
        return None

    log(
        f"  RSTE3 V3 niche cells n={len(meta_v3)} "
        f"(Tanycyte={ (meta_v3['subclass_label']=='Tanycyte NN').sum() }, "
        f"Ependymal={ (meta_v3['subclass_label']=='Ependymal NN').sum() }, "
        f"Astroependymal={ (meta_v3['subclass_label']=='Astroependymal NN').sum() })"
    )

    from bil_rstE3_fetch import read_rstE3_marker_expression

    expr_markers = read_rstE3_marker_expression(RSTE3_CURATED_MARKERS, log=log)
    present_markers = {str(g).upper() for g in expr_markers.index}
    log(f"  RSTE3 Resolve panel markers present: {len(present_markers)}/{len(RSTE3_CURATED_MARKERS)}")

    meta_v3["v3_spatial_interface"] = _v3_spatial_interface_mask(meta_v3)
    meta_v3_spatial = meta_v3.loc[meta_v3["v3_spatial_interface"]].copy()
    log(
        f"  RSTE3 V3 spatial-interface mask: n={len(meta_v3_spatial)} "
        f"({100 * len(meta_v3_spatial) / max(len(meta_v3), 1):.1f}% of subclass-filtered niche cells)"
    )

    subclass_summary = _compute_age_by_subclass_summary(meta_v3)
    subclass_path = out_dir / "exploratory_allen_rstE3_age_by_subclass_summary.csv"
    subclass_summary.to_csv(subclass_path, index=False)

    spatial_summary = _compute_age_by_subclass_summary(meta_v3_spatial)
    spatial_summary["mask"] = "v3_spatial_interface"
    spatial_path = out_dir / "exploratory_allen_rstE3_spatial_mask_age_by_subclass_summary.csv"
    spatial_summary.to_csv(spatial_path, index=False)

    s22_logfc = _load_s22_marker_logfc(out_dir)
    age_stats = _mannwhitney_gene_tests(meta_v3, expr_markers, s22_logfc)
    stats_path = out_dir / "exploratory_allen_rstE3_cell_level_age_stats.csv"
    age_stats.to_csv(stats_path, index=False)

    age_stats_spatial = _mannwhitney_gene_tests(meta_v3_spatial, expr_markers, s22_logfc)
    age_stats_spatial["spatial_mask"] = "v3_spatial_interface"
    spatial_stats_path = out_dir / "exploratory_allen_rstE3_spatial_mask_cell_level_age_stats.csv"
    age_stats_spatial.to_csv(spatial_stats_path, index=False)

    sig = age_stats[age_stats.get("fdr_bh", pd.Series(dtype=float)) <= 0.05] if "fdr_bh" in age_stats.columns else age_stats.iloc[0:0]
    n_sig_nominal = int((age_stats["mannwhitney_p"] < 0.05).sum()) if not age_stats.empty else 0
    n_sig_fdr = len(sig)

    concordant = age_stats["direction_concordant_with_s22"] == True  # noqa: E712
    n_concordant = int(concordant.sum()) if concordant.any() else 0
    n_with_s22 = int(age_stats["gse188646_s22_logFC"].notna().sum()) if not age_stats.empty else 0

    sig_spatial = (
        age_stats_spatial[age_stats_spatial.get("fdr_bh", pd.Series(dtype=float)) <= 0.05]
        if "fdr_bh" in age_stats_spatial.columns
        else age_stats_spatial.iloc[0:0]
    )
    n_sig_fdr_spatial = len(sig_spatial)

    return {
        "spatial_data_downloaded": True,
        "n_rstE3_cells_total": int(len(meta)),
        "n_rstE3_v3_niche_cells": int(len(meta_v3)),
        "n_rstE3_v3_spatial_interface_cells": int(len(meta_v3_spatial)),
        "spatial_mask_definition": (
            "Per-slide y-band between ependymal Q75 and tanycyte Q25 (ventricular interface; "
            "section-local RSTE3 coordinates)"
        ),
        "n_markers_in_resolve_panel": len(present_markers),
        "markers_in_panel": sorted(present_markers),
        "n_age_gene_tests": int(len(age_stats)),
        "n_sig_nominal_p05": n_sig_nominal,
        "n_sig_fdr05": n_sig_fdr,
        "n_spatial_mask_age_gene_tests": int(len(age_stats_spatial)),
        "n_spatial_mask_sig_fdr05": n_sig_fdr_spatial,
        "n_direction_concordant_with_gse188646_s22": n_concordant,
        "n_tests_with_s22_reference": n_with_s22,
        "top_age_effects_fdr": _records_json_safe(sig.head(10)) if not sig.empty else [],
        "top_spatial_mask_age_effects_fdr": _records_json_safe(sig_spatial.head(10)) if not sig_spatial.empty else [],
        "subclass_summary_path": subclass_path.name,
        "cell_level_stats_path": stats_path.name,
        "spatial_mask_summary_path": spatial_path.name,
        "spatial_mask_stats_path": spatial_stats_path.name,
    }


def run_allen_aging_spatial_validation(
    out_dir: Path,
    log,
    *,
    niche_strata: list[str] | None = None,
) -> None:
    log("\n=== Jin RSTE3 BIL spatial aging validation (V3 tanycytes/ependymal) ===")
    if niche_strata is None:
        niche_strata = ["6", "7", "15", "22"]

    access = _try_bil_access(log)
    access_path = out_dir / "exploratory_allen_rstE3_access_status.json"
    access_path.write_text(json.dumps(access, indent=2), encoding="utf-8")

    cell_stats = _compute_rstE3_cell_level_stats(out_dir, log)

    union_path = out_dir / "mirna_target_union_genes.csv"
    union_genes: set[str] = set()
    if union_path.is_file():
        u = pd.read_csv(union_path)
        col = "gene" if "gene" in u.columns else u.columns[0]
        union_genes = {str(g).strip().upper() for g in u[col].dropna()}

    niche_de = _load_niche_de_genes(out_dir, niche_strata)
    markers = RSTE3_CURATED_MARKERS

    rows = []
    for label, gene_set in [
        ("RSTE3_curated_V3_markers", markers),
        ("GSE188646_niche_strata_DE", niche_de),
        ("miRNA_target_union", union_genes),
    ]:
        overlap_m = gene_set & markers
        overlap_n = gene_set & niche_de
        rows.append(
            {
                "set_a": label,
                "n_genes": len(gene_set),
                "n_overlap_rstE3_markers": len(overlap_m),
                "overlap_rstE3_markers": ";".join(sorted(overlap_m)),
                "n_overlap_niche_de": len(overlap_n),
            }
        )

    universe = markers | niche_de | union_genes
    a, b = markers & universe, niche_de & universe
    tab = [[len(a & b), len(a - b)], [len(b - a), len(universe - a - b)]]
    try:
        or_, fp = stats.fisher_exact(tab, alternative="two-sided")
        fisher_or = float(or_) if np.isfinite(or_) else None
    except Exception:
        fisher_or, fp = None, float("nan")

    overlap_df = pd.DataFrame(rows)
    overlap_path = out_dir / "exploratory_allen_aging_spatial_rstE3_gene_overlap.csv"
    overlap_df.to_csv(overlap_path, index=False)

    s22_path = out_dir / "gse188646_strata" / "stratum_22_young_vs_aged_deg.csv"
    marker_detail = []
    if s22_path.is_file():
        deg = pd.read_csv(s22_path)
        deg["gene_u"] = deg["gene"].astype(str).str.upper()
        for m in sorted(markers):
            row = deg[deg["gene_u"] == m]
            if row.empty:
                marker_detail.append({"gene": m, "in_stratum22_de_table": False})
            else:
                r0 = row.iloc[0]
                marker_detail.append(
                    {
                        "gene": m,
                        "in_stratum22_de_table": True,
                        "logFC": float(r0.get("logFC", np.nan)),
                        "FDR": float(r0.get("FDR", r0.get("padj", np.nan))),
                    }
                )
    md_path = out_dir / "exploratory_allen_aging_spatial_rstE3_markers_in_s22.csv"
    if marker_detail:
        pd.DataFrame(marker_detail).to_csv(md_path, index=False)

    if cell_stats:
        interpretation = (
            "RSTE3 provides in situ Resolve spatial transcriptomics at the hypothalamic third-ventricle "
            "niche (Jin 2025 Fig. 4). Cell-level Mann-Whitney tests quantify aged (18 mo) vs young "
            "(2 mo) expression within tanycyte, ependymal, and astroependymal subclasses."
        )
        caveat = (
            "Resolve panel covers ~100 genes; not all curated V3 markers are measured. "
            "Spatial coordinates are section-local — subclass + HY ROI defines V3 proxy."
        )
    else:
        interpretation = (
            "RSTE3 is the only public in situ aging spatial dataset explicitly profiling V3 tanycytes "
            "and ependymal cells (Jin 2025 Fig. 4). Cell-level statistics require cached BIL tables."
        )
        caveat = (
            "If BIL download fails, curated marker overlap with GSE188646 is a conservative fallback only."
        )

    summ = {
        "reference_spatial": f"Jin et al. 2025 RSTE3 Resolve Molecular Cartography; BIL {BIL_DOI}",
        "bil_access": access,
        "spatial_data_downloaded": bool(cell_stats),
        "rstE3_cell_level": cell_stats or {},
        "fallback_analysis": "Curated RSTE3 V3 marker panel overlap with GSE188646 niche DE and miRNA union",
        "n_rstE3_curated_markers": len(markers),
        "n_overlap_markers_niche_de": len(markers & niche_de),
        "n_overlap_markers_mirna_union": len(markers & union_genes),
        "fisher_rstE3_markers_vs_niche_de_p": float(fp),
        "fisher_rstE3_markers_vs_niche_de_or": fisher_or,
        "interpretation": interpretation,
        "caveat": caveat,
    }
    json_path = out_dir / "exploratory_allen_aging_spatial_validation_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {access_path.name} (BIL accessible={access.get('accessible')})")
    if cell_stats:
        log(
            f"Wrote {cell_stats['subclass_summary_path']}, {cell_stats['cell_level_stats_path']}; "
            f"V3 cells={cell_stats['n_rstE3_v3_niche_cells']}, FDR<0.05 tests={cell_stats['n_sig_fdr05']}"
        )
    log(f"Wrote {overlap_path.name}; RSTE3 markers x niche DE n={len(markers & niche_de)}")
    log(f"Wrote {json_path.name}")
