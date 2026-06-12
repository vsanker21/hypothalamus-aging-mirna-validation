"""
External validation: Jin et al. 2025 Allen ABC aging scRNA (lightweight metadata only).

Uses Zeng-Aging-Mouse-WMB-taxonomy aging_degenes.csv (no expression matrix download).
Compares tanycyte/ependymal age-DE genes with:
  - htNSC miRNA target union (mirna_target_union_genes.csv)
  - MSigDB Hallmark ORA gene lists (enrichr_hallmark_ora_mirtarbase_union.csv)
  - GSE188646 third-ventricle niche stratum DE genes (clusters 6/7/15/22)

Outputs:
  exploratory_allen_aging_scrna_tanycyte_ependymal_agede.csv
  exploratory_allen_aging_scrna_hallmark_overlap.csv
  exploratory_allen_aging_scrna_gse188646_niche_gene_overlap.csv
  exploratory_allen_aging_scrna_validation_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from abc_atlas_fetch import read_abc_csv

V3_NICHE_SUPERtype_RE = r"Tanycyte|Ependymal|Hypendymal"
V3_NICHE_SUBCLASS_RE = r"Tanycyte|Ependymal|Hypendymal|Astro-Epen"
JIN_TANYCYTE_CLUSTER_IDS = {792, 793, 794, 795, 796, 797}


def _fisher_overlap(a_set: set[str], b_set: set[str], universe: set[str]) -> dict:
    a, b = a_set & universe, b_set & universe
    u = len(universe)
    if u == 0:
        return {"n_overlap": 0, "odds_ratio": None, "p_value": float("nan")}
    overlap = a & b
    # 2x2: in A & in B, in A & not B, not A & in B, neither
    tab = [
        [len(overlap), len(a - b)],
        [len(b - a), len(universe - a - b)],
    ]
    try:
        or_, p = stats.fisher_exact(tab, alternative="two-sided")
        or_v = float(or_) if np.isfinite(or_) else None
    except Exception:
        or_v, p = None, float("nan")
    return {
        "n_a": len(a),
        "n_b": len(b),
        "n_overlap": len(overlap),
        "overlap_genes": sorted(overlap),
        "odds_ratio": or_v,
        "p_value": float(p),
    }


def _load_union_genes(out_dir: Path) -> set[str]:
    p = out_dir / "mirna_target_union_genes.csv"
    if not p.is_file():
        return set()
    df = pd.read_csv(p)
    col = "gene" if "gene" in df.columns else df.columns[0]
    return {str(g).strip().upper() for g in df[col].dropna()}


def _load_hallmark_gene_sets(out_dir: Path) -> dict[str, set[str]]:
    p = out_dir / "enrichr_hallmark_ora_mirtarbase_union.csv"
    if not p.is_file():
        return {}
    df = pd.read_csv(p)
    out: dict[str, set[str]] = {}
    for _, row in df.iterrows():
        term = str(row.get("Term", "")).strip()
        genes_raw = str(row.get("Genes", ""))
        if not term or genes_raw == "nan":
            continue
        genes = {g.strip().upper() for g in genes_raw.replace(";", " ").split() if g.strip()}
        if genes:
            out[term] = genes
    return out


def _deg_gene_set(deg: pd.DataFrame, padj_max: float = 0.05, abs_lfc_min: float = 0.25) -> set[str]:
    genes: set[str] = set()
    padj_col = next(
        (c for c in ("FDR", "padj", "p_val_adj", "adj.P.Val", "adjusted_pvalue") if c in deg.columns),
        None,
    )
    for _, row in deg.iterrows():
        g = str(row.get("gene", "")).strip().upper()
        if not g:
            continue
        lfc = float(row.get("logFC", row.get("log2FoldChange", np.nan)))
        if padj_col:
            padj = float(row[padj_col])
            if np.isfinite(padj) and padj <= padj_max and np.isfinite(lfc) and abs(lfc) >= abs_lfc_min:
                genes.add(g)
        elif np.isfinite(lfc) and abs(lfc) >= abs_lfc_min:
            genes.add(g)
    return genes


def _load_gse188646_niche_de_genes(out_dir: Path, niche_strata: list[str]) -> dict[str, set[str]]:
    strata_dir = out_dir / "gse188646_strata"
    result: dict[str, set[str]] = {}
    for sid in niche_strata:
        deg_path = strata_dir / f"stratum_{sid}_young_vs_aged_deg.csv"
        if not deg_path.is_file():
            continue
        result[sid] = _deg_gene_set(pd.read_csv(deg_path))
    return result


def _filter_v3_astro_epen(degenes: pd.DataFrame) -> pd.DataFrame:
    d = degenes.copy()
    d["gene_symbol_u"] = d["gene_symbol"].astype(str).str.strip().str.upper()
    name = d["grouping_name"].astype(str)
    gtype = d["grouping_type"].astype(str).str.lower()
    mask_super = name.str.contains(V3_NICHE_SUPERtype_RE, case=False, na=False, regex=True)
    mask_sub = name.str.contains(V3_NICHE_SUBCLASS_RE, case=False, na=False, regex=True)
    # Jin tanycyte clusters 792-797 at cluster level
    cluster_alias = pd.to_numeric(d.get("grouping_label", pd.Series(dtype=float)), errors="coerce")
    mask_cluster = gtype.eq("cluster") & cluster_alias.isin(list(JIN_TANYCYTE_CLUSTER_IDS))
    d["v3_niche_celltype"] = mask_super | mask_sub | mask_cluster
    return d[d["v3_niche_celltype"]].copy()


def run_allen_aging_scrna_validation(
    out_dir: Path,
    log,
    *,
    niche_strata: list[str] | None = None,
) -> None:
    log("\n=== Allen Jin 2025 aging scRNA external validation (lightweight ABC metadata) ===")
    if niche_strata is None:
        niche_strata = ["6", "7", "15", "22"]

    degenes = read_abc_csv("aging_degenes", log=log)
    if degenes is None or degenes.empty:
        log("Allen aging scRNA validation: aging_degenes.csv unavailable; skipped.")
        return

    v3_de = _filter_v3_astro_epen(degenes)
    if v3_de.empty:
        log("Allen aging scRNA validation: no tanycyte/ependymal age-DE rows matched; skipped.")
        return

    v3_de_out = out_dir / "exploratory_allen_aging_scrna_tanycyte_ependymal_agede.csv"
    v3_de.to_csv(v3_de_out, index=False)

    # Unique age-DE genes across V3 niche types (any significant row)
    v3_gene_set = set(v3_de["gene_symbol_u"].dropna())
    union_set = _load_union_genes(out_dir)
    hallmark_sets = _load_hallmark_gene_sets(out_dir)
    niche_de = _load_gse188646_niche_de_genes(out_dir, niche_strata)

    # Universe: all genes in aging_degenes table
    universe = {str(g).strip().upper() for g in degenes["gene_symbol"].dropna()}

    union_overlap = _fisher_overlap(union_set, v3_gene_set, universe)

    # Hallmark overlaps
    hm_rows = []
    for term, genes in hallmark_sets.items():
        ov = _fisher_overlap(genes, v3_gene_set, universe)
        hm_rows.append(
            {
                "hallmark_term": term,
                "n_hallmark_genes": ov["n_a"],
                "n_jin_v3_agede": ov["n_b"],
                "n_overlap": ov["n_overlap"],
                "odds_ratio": ov["odds_ratio"],
                "fisher_p": ov["p_value"],
                "overlap_genes": ";".join(ov["overlap_genes"][:50]),
            }
        )
    hm_df = pd.DataFrame(hm_rows).sort_values("fisher_p")
    hm_path = out_dir / "exploratory_allen_aging_scrna_hallmark_overlap.csv"
    hm_df.to_csv(hm_path, index=False)

    # GSE188646 niche stratum DE overlap
    niche_rows = []
    all_niche_de: set[str] = set()
    for sid, genes in niche_de.items():
        all_niche_de |= genes
        ov_u = _fisher_overlap(genes, v3_gene_set, universe)
        ov_union = _fisher_overlap(genes, union_set, universe)
        niche_rows.append(
            {
                "gse188646_stratum": sid,
                "n_stratum_de_genes": len(genes),
                "n_overlap_jin_v3_agede": ov_u["n_overlap"],
                "fisher_p_jin_v3": ov_u["p_value"],
                "n_overlap_mirna_union": ov_union["n_overlap"],
                "fisher_p_union": ov_union["p_value"],
            }
        )
    niche_df = pd.DataFrame(niche_rows)
    niche_path = out_dir / "exploratory_allen_aging_scrna_gse188646_niche_gene_overlap.csv"
    niche_df.to_csv(niche_path, index=False)

    # Direction concordance for genes in both Jin tanycyte DE and GSE188646 stratum 22 DE
    direction_rows = []
    s22 = niche_de.get("22", set())
    shared_dir = s22 & v3_gene_set
    if shared_dir and (out_dir / "gse188646_strata" / "stratum_22_young_vs_aged_deg.csv").is_file():
        deg22 = pd.read_csv(out_dir / "gse188646_strata" / "stratum_22_young_vs_aged_deg.csv")
        deg22["gene_u"] = deg22["gene"].astype(str).str.upper()
        jin_tany = v3_de[v3_de["grouping_name"].astype(str).str.contains("Tanycyte", case=False, na=False)]
        jin_by_gene = jin_tany.groupby("gene_symbol_u")["age_effect_size"].median()
        for g in sorted(shared_dir):
            if g not in jin_by_gene.index:
                continue
            row = deg22[deg22["gene_u"] == g]
            if row.empty:
                continue
            gse_lfc = float(row.iloc[0]["logFC"])
            jin_eff = float(jin_by_gene[g])
            direction_rows.append(
                {
                    "gene": g,
                    "gse188646_stratum22_logFC": gse_lfc,
                    "jin_tanycyte_age_effect_size": jin_eff,
                    "same_direction": bool(np.sign(gse_lfc) == np.sign(jin_eff)),
                }
            )
    dir_df = pd.DataFrame(direction_rows)
    if not dir_df.empty:
        dir_path = out_dir / "exploratory_allen_aging_scrna_direction_concordance_s22.csv"
        dir_df.to_csv(dir_path, index=False)
        n_conc = int(dir_df["same_direction"].sum())
        frac_conc = n_conc / len(dir_df) if len(dir_df) else float("nan")
    else:
        n_conc, frac_conc = 0, float("nan")

    # Count age-DE by grouping level for V3 types
    de_counts = (
        v3_de.groupby(["grouping_type", "grouping_name"])["gene_symbol_u"]
        .nunique()
        .reset_index(name="n_age_de_genes")
        .sort_values("n_age_de_genes", ascending=False)
    )

    summ = {
        "reference": "Jin et al. 2025 Nature; Allen ABC Zeng-Aging-Mouse-WMB-taxonomy aging_degenes.csv",
        "data_access": "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/aging_degenes.csv",
        "n_total_agede_rows": int(len(degenes)),
        "n_v3_niche_agede_rows": int(len(v3_de)),
        "n_unique_v3_niche_agede_genes": len(v3_gene_set),
        "n_mirna_union_genes": len(union_set),
        "mirna_union_vs_jin_v3_agede_overlap": union_overlap["n_overlap"],
        "mirna_union_vs_jin_v3_agede_fisher_p": union_overlap["p_value"],
        "mirna_union_vs_jin_v3_agede_odds_ratio": union_overlap["odds_ratio"],
        "top_hallmark_jin_v3_overlaps": hm_df.head(5)[
            ["hallmark_term", "n_overlap", "fisher_p"]
        ].to_dict(orient="records"),
        "gse188646_niche_strata": niche_strata,
        "gse188646_niche_de_overlap_rows": niche_rows,
        "n_shared_genes_s22_jin_tanycyte_direction_test": len(direction_rows),
        "direction_concordance_fraction_s22": frac_conc,
        "top_v3_celltypes_by_n_agede": de_counts.head(10).to_dict(orient="records"),
        "interpretation": (
            "Independent brain-wide aging atlas (Jin 2025) reports extensive age-DE in tanycytes "
            "and ependymal cells at the third ventricle — supporting anatomical/biological plausibility "
            "of the GSE188646 niche panel. Overlap with the miRNA target union tests whether literature-prior "
            "targets co-occur with V3 aging signatures (enrichment ≠ causal miRNA regulation). "
            "Hallmark overlap links V3 aging genes to inflammatory/stress axes already seen in ORA."
        ),
        "caveat": (
            "Jin cohort: 2 vs 18 mo, both sexes, broad brain dissection including hypothalamus — "
            "not female-only snRNA like GSE188646. age_effect_size sign/magnitude differs from logFC scale."
        ),
    }
    json_path = out_dir / "exploratory_allen_aging_scrna_validation_summary.json"
    json_path.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    log(f"Wrote {v3_de_out.name} ({len(v3_de)} V3 niche age-DE rows).")
    log(f"Wrote {hm_path.name}; union∩Jin-V3 overlap n={union_overlap['n_overlap']}, Fisher P={union_overlap['p_value']:.3g}".replace("∩", " intersect "))
    log(f"Wrote {json_path.name}")
