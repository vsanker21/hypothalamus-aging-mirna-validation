"""
Science Advances–oriented computational suite: hypothalamic niche / NSC-like biology,
aging pseudobulk, miRNA-target structure, and LIFU-relevant pathway priors.

All outputs use prefix exploratory_sa_nsc_lifu_* (exploratory; not confirmatory).

Gene sets are small literature priors (mouse symbols); intersected with assay genes before use.
Cite original papers when elevating beyond hypothesis generation.
"""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

from config import OUTPUT_DIR
from mirna_target_union import ensure_mirtarbase_gmt, mmu_gmt_gene_universe
from overlap_fisher import fisher_overlap

# Mouse symbols (mixed sources: neural stem / radial glia–like, tanycyte–ependymal, niche glia,
# mechanosensitive / calcium / neuroinflammatory priors relevant to ultrasound biophysics narratives).
# Intersection with data removes absent symbols.
CURATED_GENE_SETS: dict[str, list[str]] = {
    "NSC_like_core": [
        "Sox2",
        "Nes",
        "Hes1",
        "Hes5",
        "Ascl1",
        "Prom1",
        "Fabp7",
        "Gfap",
        "Slc1a3",
        "Aqp4",
    ],
    "Radial_glia_like": ["Vim", "Rbpms", "Pax6", "Hes1", "Hes5", "Fabp7", "Gfap", "Slc1a3"],
    "Tanycyte_ependymal_like": ["Rax", "Dcdc2a", "Col23a1", "Ccdc153", "Foxj1", "Pkd2l1", "Dnah5"],
    "OPC_oligolineage": ["Pdgfra", "Cspg4", "Olig2", "Sox10", "Gpr17", "Bcas1"],
    "Microglia_homeostasis": ["Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Hexb", "C1qa", "C1qb"],
    "LIFU_mechanosensory_calcium": ["Piezo1", "Piezo2", "Trpv4", "Trpv1", "Itpr2", "P2rx7", "Calm1", "Camk2a"],
    "Neuroinflammatory_signaling": ["Tnf", "Il1b", "Nfkb1", "Nlrp3", "Tlr4", "Rela", "Stat3", "Il6"],
}


def _norm_sym(s: str) -> str:
    return str(s).strip().upper()


def _spearman_xy(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    r = stats.spearmanr(x, y, nan_policy="omit")
    stat = float(getattr(r, "statistic", getattr(r, "correlation", np.nan)))
    pv = float(getattr(r, "pvalue", np.nan))
    return stat, pv


def _write_curated_gmt(path: Path, sets: dict[str, list[str]], valid: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for name, genes in sets.items():
        g2 = sorted({_norm_sym(g) for g in genes if _norm_sym(g) in valid})
        if len(g2) < 5:
            continue
        safe = re.sub(r"[^\w]+", "_", name).strip("_")
        lines.append(safe + "\tna\t" + "\t".join(g2))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _pseudobulk_module_scores(
    mtx_path: Path,
    row_path: Path,
    col_path: Path,
    meta_path: Path,
    sets: dict[str, list[str]],
) -> pd.DataFrame | None:
    if not all(p.is_file() for p in (mtx_path, row_path, col_path, meta_path)):
        return None
    X = mmread(mtx_path).tocsr()
    rows = pd.read_csv(row_path).iloc[:, 0].astype(str).str.strip().str.upper().tolist()
    cols = pd.read_csv(col_path).iloc[:, 0].astype(str).str.strip().tolist()
    meta = pd.read_csv(meta_path)
    if "orig.ident" not in meta.columns or "age_bin" not in meta.columns:
        return None
    col_idx = {c: i for i, c in enumerate(cols)}
    idx_y = [col_idx[c] for c in cols if c in col_idx and meta.set_index("orig.ident").loc[c, "age_bin"] == "Young"]
    idx_a = [col_idx[c] for c in cols if c in col_idx and meta.set_index("orig.ident").loc[c, "age_bin"] == "Aged"]
    if not idx_y or not idx_a:
        return None
    # columns are samples; normalize to log1p CPM per column
    Xc = X.tocsc()
    dense_cols = []
    for j in range(Xc.shape[1]):
        col = Xc.getcol(j).toarray().ravel().astype(float)
        s = col.sum()
        if s <= 0:
            dense_cols.append(np.zeros_like(col))
        else:
            dense_cols.append(np.log1p(col / s * 1e6))
    M = np.column_stack(dense_cols)  # genes x samples
    row_to_i = {g: i for i, g in enumerate(rows)}
    rows_out = []
    for set_name, genes in sets.items():
        gi = [row_to_i[_norm_sym(g)] for g in genes if _norm_sym(g) in row_to_i]
        if len(gi) < 3:
            continue
        v = M[gi, :].mean(axis=0)
        yv = v[idx_y]
        av = v[idx_a]
        tt = stats.ttest_ind(yv, av, equal_var=False)
        rows_out.append(
            {
                "gene_set": set_name,
                "n_genes_used": len(gi),
                "mean_score_young": float(np.mean(yv)),
                "mean_score_aged": float(np.mean(av)),
                "welch_t": float(tt.statistic),
                "p_two_sided": float(tt.pvalue),
            }
        )
    return pd.DataFrame(rows_out) if rows_out else None


def _fgsea_custom(
    ranked: pd.DataFrame,
    gmt_path: Path,
    outdir: Path,
) -> pd.DataFrame | None:
    try:
        import gseapy as gp
    except ImportError:
        return None
    rnk = ranked[["gene", "rank_metric"]].dropna().drop_duplicates(subset=["gene"])
    if len(rnk) < 50:
        return None
    rnk = rnk.set_index("gene")
    rnk.columns = ["rank"]
    outdir.mkdir(parents=True, exist_ok=True)
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=str(gmt_path),
        outdir=str(outdir),
        permutation_num=499,
        seed=11,
        verbose=False,
        no_plot=True,
        min_size=5,
        max_size=2000,
    )
    df = pre.res2d
    return df if df is not None and len(df) else None


def run_suite(out_dir: Path | None = None, log=print) -> None:
    out_dir = out_dir or OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    log("\n=== Science Advances framing suite: NSC/niche + aging + LIFU priors (exploratory) ===")

    deg_path = out_dir / "gse188646_young_vs_aged_deg.csv"
    if not deg_path.is_file():
        log("SA suite: missing gse188646_young_vs_aged_deg.csv; skipped.")
        return
    deg = pd.read_csv(deg_path)
    deg["gene_u"] = deg["gene"].astype(str).map(_norm_sym)
    universe = set(deg["gene_u"])

    gmt_mmu = ensure_mirtarbase_gmt()
    uni_mmu = {g.upper() for g in mmu_gmt_gene_universe(gmt_mmu)}
    union_path = out_dir / "mirna_target_union_genes.csv"
    union: set[str] = set()
    if union_path.is_file():
        udf = pd.read_csv(union_path)
        col = "gene" if "gene" in udf.columns else udf.columns[0]
        union = {str(x).strip().upper() for x in udf[col].dropna()}

    # valid symbols: in DE table (pseudobulk universe)
    valid_for_gmt = universe

    gmt_custom = out_dir / "sa_nsc_lifu_curated_sets.gmt"
    _write_curated_gmt(gmt_custom, CURATED_GENE_SETS, valid_for_gmt)

    # Rank metric for prerank GSEA
    pv = np.maximum(deg["p_val"].astype(float).values, 1e-300)
    deg["rank_metric"] = np.sign(deg["logFC"].astype(float)) * (-np.log10(pv))
    ranked = deg[["gene", "rank_metric"]].copy()
    ranked["gene"] = ranked["gene"].astype(str).str.strip().str.upper()

    gsea_dir = out_dir / "sa_nsc_lifu_gsea_prerank_curated"
    gsea_df = _fgsea_custom(ranked, gmt_custom, gsea_dir)
    if gsea_df is not None:
        gsea_df.to_csv(out_dir / "exploratory_sa_nsc_lifu_fgsea_curated_sets.csv", index=False)
        log(f"Wrote exploratory_sa_nsc_lifu_fgsea_curated_sets.csv ({len(gsea_df)} rows)")

    # Fisher: each curated set vs miRNA union inside miRTarBase universe
    fish_rows = []
    for name, genes in CURATED_GENE_SETS.items():
        ext = {_norm_sym(g) for g in genes} & uni_mmu
        if len(ext) < 5 or not union:
            continue
        res = fisher_overlap(union, ext, uni_mmu)
        res["gene_set"] = name
        fish_rows.append(res)
    if fish_rows:
        pd.DataFrame(fish_rows).to_csv(out_dir / "exploratory_sa_nsc_lifu_fisher_targets_vs_curated_sets.csv", index=False)
        log("Wrote exploratory_sa_nsc_lifu_fisher_targets_vs_curated_sets.csv")

    # Pseudobulk module young vs aged
    pb = _pseudobulk_module_scores(
        out_dir / "gse188646_pseudobulk_counts.mtx",
        out_dir / "gse188646_pseudobulk_counts_rownames.csv",
        out_dir / "gse188646_pseudobulk_counts_colnames.csv",
        out_dir / "gse188646_pseudobulk_metadata.csv",
        CURATED_GENE_SETS,
    )
    if pb is not None:
        pb.to_csv(out_dir / "exploratory_sa_nsc_lifu_pseudobulk_module_young_vs_aged.csv", index=False)
        log("Wrote exploratory_sa_nsc_lifu_pseudobulk_module_young_vs_aged.csv")

    # Cross-modal burden within niche gene sets vs rest
    bur_path = out_dir / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv"
    if bur_path.is_file():
        bdf = pd.read_csv(bur_path)
        bdf["gene_u"] = bdf["gene"].astype(str).map(_norm_sym)
        niche_union = set()
        for genes in CURATED_GENE_SETS.values():
            niche_union |= {_norm_sym(g) for g in genes}
        bdf["in_curated_niche_union"] = bdf["gene_u"].isin(niche_union)
        r_all = _spearman_xy(bdf["weighted_burden"], bdf["logFC"])
        sub = bdf[bdf["in_curated_niche_union"]]
        rest = bdf[~bdf["in_curated_niche_union"]]
        r_sub = _spearman_xy(sub["weighted_burden"], sub["logFC"]) if len(sub) > 15 else (float("nan"), float("nan"))
        r_rest = _spearman_xy(rest["weighted_burden"], rest["logFC"]) if len(rest) > 15 else (float("nan"), float("nan"))
        pd.DataFrame(
            [
                {
                    "subset": "all_genes_with_burden",
                    "n": int(len(bdf)),
                    "spearman_rho": r_all[0],
                    "spearman_p": r_all[1],
                },
                {
                    "subset": "in_curated_niche_union",
                    "n": int(len(sub)),
                    "spearman_rho": r_sub[0],
                    "spearman_p": r_sub[1],
                },
                {
                    "subset": "not_in_niche_union",
                    "n": int(len(rest)),
                    "spearman_rho": r_rest[0],
                    "spearman_p": r_rest[1],
                },
            ]
        ).to_csv(out_dir / "exploratory_sa_nsc_lifu_burden_logfc_by_niche_subset.csv", index=False)
        log("Wrote exploratory_sa_nsc_lifu_burden_logfc_by_niche_subset.csv")

    # HypoMap CELLxGENE Spearman subset: niche-like reference types
    hypo_c = out_dir / "gse188646_hypomap_mapping" / "hypomap_custom_ref_spearman.csv"
    if hypo_c.is_file():
        h = pd.read_csv(hypo_c)
        pat = re.compile(r"OPC|NSC|Stem|Astrocyte|Tanycyte|Ependymal|Radial|Neural", re.I)
        h["niche_like"] = h["ref_celltype"].astype(str).apply(lambda s: bool(pat.search(s)))
        h.to_csv(out_dir / "exploratory_sa_nsc_lifu_hypomap_ref_flags.csv", index=False)
        agg = h.groupby("niche_like")["rho"].agg(["mean", "median", "count"]).reset_index()
        agg.to_csv(out_dir / "exploratory_sa_nsc_lifu_hypomap_rho_summary_by_niche_flag.csv", index=False)
        log("Wrote exploratory_sa_nsc_lifu_hypomap_ref_flags.csv")

    # Cluster annotation flags
    lab = out_dir / "gse188646_cluster_annotation" / "cluster_putative_labels.csv"
    if lab.is_file():
        L = pd.read_csv(lab)
        pat2 = re.compile(r"OPC|Tanycyte|Astrocyte|Microglia|Ependymal", re.I)
        L["niche_glial_flag"] = (
            L["rank1_module"].astype(str).apply(lambda s: bool(pat2.search(s)))
            | L["rank2_module"].astype(str).apply(lambda s: bool(pat2.search(s)))
        )
        L.to_csv(out_dir / "exploratory_sa_nsc_lifu_cluster_niche_flags.csv", index=False)
        log("Wrote exploratory_sa_nsc_lifu_cluster_niche_flags.csv")

    summ = {
        "completed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "outputs_prefix": "exploratory_sa_nsc_lifu_*",
        "gmt_written": str(gmt_custom.resolve()) if gmt_custom.is_file() else None,
        "methodology_note": "Prerank GSEA uses signed -log10(p) * sign(logFC) on DE genes; gene index uppercased for gseapy GMT matching. Fisher uses miRTarBase mmu universe. Pseudobulk modules: mean log1p-CPM per set, Welch Young vs Aged.",
        "note": "Curated gene sets are priors for hypothesis generation; intersect with DE universe for GSEA/Fisher.",
    }
    (out_dir / "exploratory_sa_nsc_lifu_suite_summary.json").write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log("Wrote exploratory_sa_nsc_lifu_suite_summary.json")
    log("SA suite: done.")
