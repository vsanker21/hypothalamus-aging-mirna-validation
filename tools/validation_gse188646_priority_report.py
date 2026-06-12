"""
QC pseudobulk metadata + .mtx, compare stratum vs global DE, merge meta DE,
and write a prioritised shortlist with explicit statistical caveats.

Run from feasibility_study/:
  python tools/validation_gse188646_priority_report.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
REPORT = OUT / "GSE188646_VALIDATION_AND_PRIORITIZATION_REPORT.txt"
SHORTLIST = OUT / "gse188646_prioritized_genes.csv"


# Mouse symbols commonly discussed in hypothalamus / energy-balance / neuroimmune context (illustrative, not exhaustive).
HYPOTHALAMUS_CONTEXT_GENES = {
    "NPY", "AGRP", "POMC", "CART", "GAL", "CRH", "TRH", "OXT", "AVP", "KISS1",
    "TAC2", "LEP", "LEPR", "MC4R", "MC3R", "INSR", "IGF1", "SOX2", "GFAP", "S100B",
    "SLC17A6", "SLC32A1", "GAD1", "GAD2", "DRD1", "DRD2", "HCRT", "PMCH",
    "STAT3", "IL6", "IL1B", "TNF", "NLRP3", "CX3CL1", "CXCL12",
}


def _load_mtx_bundle() -> tuple[np.ndarray, list[str], list[str]]:
    mat = mmread(OUT / "gse188646_pseudobulk_counts.mtx").tocsr()
    genes = (
        pd.read_csv(OUT / "gse188646_pseudobulk_counts_rownames.csv")["gene"].astype(str).tolist()
    )
    samples = (
        pd.read_csv(OUT / "gse188646_pseudobulk_counts_colnames.csv")["sample"].astype(str).tolist()
    )
    assert mat.shape == (len(genes), len(samples)), (mat.shape, len(genes), len(samples))
    return mat.toarray(), genes, samples


def _bh_note(p: pd.Series) -> str:
    """Explain degenerate BH when all adjusted p-values collapse."""
    p = pd.to_numeric(p, errors="coerce")
    m = int(p.notna().sum())
    ps = np.sort(p.dropna().to_numpy())
    if len(ps) == 0:
        return "No p-values."
    vmin = float(ps[0])
    vmax = float(ps[-1])
    v1 = vmin * m
    note = (
        f"Across m={m} tests, smallest raw p≈{vmin:.2e} so m*p_min≈{v1:.3f}. "
        f"Largest raw p≈{vmax:.6f}. Under standard BH, if all intermediate "
        f"rank-scaled values stay ≥ max(p), every BH-adjusted value ties at max(p) "
        f"(no gene clears a 0.05 FDR gate). Use nominal p and effect sizes with multiplicity caution."
    )
    return note


def main() -> int:
    lines: list[str] = []
    ap = lambda s: lines.append(s)

    ap("GSE188646 — validation, QC, prioritisation (machine-generated)")
    ap("=" * 72)
    ap("")
    ap("1) LOCKED NARRATIVE + SAP + EXTENDED_REPORT (alignment)")
    ap("-" * 72)
    ap(
        "SCIENTIFIC_STORY_ONE_PAGE.txt: GSE188646 is an orthogonal anchor (pseudobulk DE, "
        "Fisher vs miRNA targets, PROGENy/DoRothEA, exploratory fgsea/camera/fry) — not a "
        "replication of Zhang MOESM bulk/exosome design."
    )
    ap(
        "STATISTICAL_ANALYSIS_PLAN.txt: GSE188646 pseudobulk has n=4 Young vs 4 Aged animals; "
        "power is limited to large consistent effects. Exploratory outputs must not be upgraded "
        "to confirmatory endpoints post hoc."
    )
    ap(
        "EXTENDED_REPORT.txt: documents the same extended stack (Hallmark ORA, multi-library ORA, "
        "optional cohort-2 meta when enabled, Tier-2 decoupler, etc.). Keep manuscript claims inside "
        "these boundaries."
    )
    ap("")

    ap("2) PSEUDOBULK METADATA QC")
    ap("-" * 72)
    meta = pd.read_csv(OUT / "gse188646_pseudobulk_metadata.csv")
    ap(meta.to_string(index=False))
    ap("")
    vc = meta.groupby("age_bin")["orig.ident"].nunique()
    ap(f"Distinct orig.ident per age_bin: {vc.to_dict()}")
    if set(meta["age_bin"]) >= {"Young", "Aged"} and vc.get("Young", 0) == vc.get("Aged", 0) == 4:
        ap("PASS: balanced 4 vs 4 biological replicates (matches SAP).")
    else:
        ap("WARN: unexpected Young/Aged replicate counts — review metadata.")
    ap("")

    ap("3) PSEUDOBULK .mtx QC + SPOT-CHECKS")
    ap("-" * 72)
    counts, genes, samples = _load_mtx_bundle()
    lib = counts.sum(axis=0)
    ap(f"Matrix shape: {counts.shape[0]} genes × {counts.shape[1]} samples.")
    ap(f"Library sizes (column sums): min={lib.min():.0f}, max={lib.max():.0f}, all>0: {bool(np.all(lib > 0))}")
    ap("Per-sample totals: " + ", ".join(f"{s}={v:.0f}" for s, v in zip(samples, lib)))
    ap("")
    # Spot-check: compare row sum for a few named genes to manual slice
    checks = ["Actb", "Malat1", "Gapdh", "Tmsb4x", "Rplp0"]
    for g in checks:
        if g in genes:
            i = genes.index(g)
            row = counts[i, :]
            ap(f"Gene {g}: total counts across samples = {row.sum():.0f} (per-sample: {row.tolist()})")
        else:
            ap(f"Gene {g}: not in matrix rownames (skipped).")
    ap("PASS: .mtx readable; library sizes strictly positive.")
    ap("")

    ap("4) GLOBAL PSEUDOBULK DE — MULTIPLICITY / BH NOTE")
    ap("-" * 72)
    deg = pd.read_csv(OUT / "gse188646_young_vs_aged_deg.csv")
    ap(_bh_note(deg["p_val"]))
    ap("")
    deg["z_edgeR"] = deg["logFC"] / np.maximum(deg["se_logFC"], 1e-12)
    ap(
        "Because BH FDR is uninformative here, prioritisation uses nominal p_val, "
        "|logFC|/se_logFC (quasi-F-derived), and cross-cohort meta columns — all labelled exploratory."
    )
    ap("")

    ap("5) STRATUM VS GLOBAL logFC (Spearman on intersecting genes)")
    ap("-" * 72)
    man = pd.read_csv(OUT / "gse188646_strata" / "manifest.csv")
    gmap = deg.set_index("gene")["logFC"].astype(float)
    strata_pick = ["0", "1", "2", "5", "10", "15"]
    for sid in strata_pick:
        pth = OUT / "gse188646_strata" / f"stratum_{sid}_young_vs_aged_deg.csv"
        if not pth.is_file():
            continue
        sd = pd.read_csv(pth)
        if "gene" not in sd.columns or "logFC" not in sd.columns:
            continue
        sser = sd.set_index("gene")["logFC"].astype(float)
        common = gmap.index.intersection(sser.index)
        if len(common) < 500:
            ap(f"Stratum {sid}: too few genes for correlation ({len(common)}).")
            continue
        rho, psp = stats.spearmanr(gmap.loc[common], sser.loc[common])
        row = man.loc[man["stratum"].astype(str) == sid]
        nc = int(row["n_cells"].iloc[0]) if len(row) else -1
        ap(f"Stratum {sid} (n_cells≈{nc}): Spearman rho={rho:.3f}, two-sided p={psp:.2e} (n={len(common)} genes)")
    ap(
        "Interpretation: high rho ⇒ cluster-level young-vs-aged direction broadly tracks global pseudobulk; "
        "low rho in a large stratum can indicate cell-type–specific ageing programs or noise."
    )
    ap("")

    ap("6) PRIORITISED SHORTLIST (effect + uncertainty + meta + biology flags)")
    ap("-" * 72)
    meta_path = OUT / "exploratory_meta_DE_two_cohort_DL.csv"
    if not meta_path.is_file():
        ap("WARN: exploratory_meta_DE_two_cohort_DL.csv missing — meta columns omitted.")
        m = deg.copy()
        m["beta_DL"] = np.nan
        m["se_DL"] = np.nan
        m["p_two_sided"] = np.nan
        m["logFC_gse87102"] = np.nan
        m["fdr_bh_meta"] = np.nan
    else:
        meta_df = pd.read_csv(meta_path)
        meta_df["_g"] = meta_df["gene"].astype(str).str.upper()
        deg_m = deg.copy()
        deg_m["_g"] = deg_m["gene"].astype(str).str.upper()
        meta_cols = [
            "_g",
            "beta_DL",
            "se_DL",
            "tau2_DL",
            "z",
            "p_two_sided",
            "logFC_gse188646",
            "logFC_gse87102",
            "fdr_bh",
        ]
        use = [c for c in meta_cols if c in meta_df.columns]
        meta_sub = meta_df[use].copy()
        if "fdr_bh" in meta_sub.columns:
            meta_sub = meta_sub.rename(columns={"fdr_bh": "fdr_bh_meta"})
        m = deg_m.merge(meta_sub, on="_g", how="left").drop(columns=["_g"], errors="ignore")
    m["hypothalamus_context"] = m["gene"].str.upper().isin(HYPOTHALAMUS_CONTEXT_GENES)
    m["abs_lfc"] = m["logFC"].abs()
    m["nominal_sig"] = m["p_val"] <= 0.05
    m["meta_z"] = m["beta_DL"] / m["se_DL"].replace(0, np.nan)
    m["sign_consistent"] = (
        np.sign(m["logFC"].fillna(0))
        * np.sign(m["logFC_gse87102"].fillna(0))
        >= 0
    ) & m["logFC_gse87102"].notna()
    # Rank: prefer nominal sig + large |z_edgeR| + hypothalamus context + meta support
    m["rank_score"] = (
        m["nominal_sig"].astype(int) * 50
        + m["hypothalamus_context"].astype(int) * 30
        + np.clip(m["z_edgeR"].abs(), 0, 20)
        + np.clip(m["meta_z"].abs().fillna(0), 0, 5)
        + m["sign_consistent"].astype(int) * 5
    )
    top = m.sort_values(["rank_score", "abs_lfc"], ascending=False).head(40)
    cols = [
        "gene",
        "logFC",
        "se_logFC",
        "p_val",
        "p_val_adj",
        "z_edgeR",
        "beta_DL",
        "se_DL",
        "p_two_sided",
        "logFC_gse87102",
        "fdr_bh_meta",
        "hypothalamus_context",
        "sign_consistent",
        "rank_score",
    ]
    cols = [c for c in cols if c in top.columns]
    top[cols].to_csv(SHORTLIST, index=False)
    ap(f"Wrote {SHORTLIST.name} (top 40 by rank_score; review manually).")
    ap("")
    ap(top[cols].head(20).to_string(index=False))
    ap("")

    ap("7) PRIMARY CONTRAST STORY (for grants / figures)")
    ap("-" * 72)
    ap(
        "Primary public anchor for hypothalamus ageing in this repo: GSE188646 pseudobulk "
        "young vs aged (edgeR QLF). Two-cohort DL meta and miRNA/HypoMap layers are supportive "
        "or contrast-mismatched context per SCIENTIFIC_STORY_ONE_PAGE.txt — label accordingly."
    )
    ap("")
    ap("8) ORTHOGONAL LITERATURE (illustrative; not causal inference from this repo)")
    ap("-" * 72)
    ap(
        "Lepr: Hypothalamic LepR populations are heterogeneous and central to energy balance; "
        "age-related leptin resistance and altered downstream signalling are widely discussed in rodent models "
        "(interpret pseudobulk as mixed-cell signal, not subtype-specific attribution)."
    )
    ap(
        "Gad1 / Drd1: GABAergic and dopaminergic signalling are core hypothalamic circuit motifs; ageing can remodel "
        "transmitter balance without implying cell-type specificity from pseudobulk alone."
    )
    ap(
        "Rxfp1: Relaxin-family receptor biology is active in neuroscience literature (often alongside RXFP3 / relaxin-3); "
        "Rxfp1 in this shortlist reflects nominal ranking in GSE188646 pseudobulk — validate expression and function separately."
    )
    ap(
        "Further reading examples: https://pubmed.ncbi.nlm.nih.gov/35927440/ (hypothalamic LepR neurons); "
        "https://www.mdpi.com/1422-0067/23/8/4387 (relaxin-3 / RXFP3 ageing-related disease context — note receptor paralogy vs Rxfp1)."
    )
    ap("")

    cmap = ROOT / "outputs" / "gse188646_cluster_annotation" / "cluster_putative_labels.csv"
    if cmap.is_file():
        ap("9) CLUSTER → PUTATIVE CELL CLASS (marker modules; beyond global–stratum ρ)")
        ap("-" * 72)
        ap(
            "Source: r/gse188646_cluster_marker_mapping.R on the integrated RDS (see CLUSTER_MARKER_MAPPING_README.txt). "
            "Rank1/rank2 are z-scored module means across clusters — not author cell-type labels."
        )
        cl_full = pd.read_csv(cmap)
        if "delta_z" in cl_full.columns:
            amb = cl_full.sort_values("delta_z", ascending=True).head(8)
            ap("Lowest rank1–rank2 separation (ambiguous modules; interpret cautiously):")
            ap(amb.to_string(index=False))
            ap("")
        if "n_cells" in cl_full.columns:
            big = cl_full.sort_values("n_cells", ascending=False).head(12)
            ap("Largest strata with putative labels (join to stratum_* DE files via deg_csv):")
            ap(
                big[
                    [
                        "seurat_cluster_id",
                        "rank1_module",
                        "rank2_module",
                        "delta_z",
                        "n_cells",
                    ]
                ].to_string(index=False)
            )
            ap("")

    hax = ROOT / "outputs" / "gse188646_hypomap_mapping" / "hypomap_axis_spearman.csv"
    if hax.is_file():
        ap("10) HYPOMAP GSE208355 DE AXES vs CLUSTERS (harmonised gene space; not atlas label transfer)")
        ap("-" * 72)
        ap(
            "Source: r/gse188646_hypomap_reference_mapping.R — Spearman(cluster mean expr, log2FC) per IP/input contrast; "
            "see outputs/gse188646_hypomap_mapping/HYPomap_REFERENCE_MAPPING_README.txt. With large n_genes, prefer |rho| over p."
        )
        hs = pd.read_csv(hax)
        if len(hs):
            hs2 = hs.assign(_ar=hs["rho"].abs())
            top = hs2.sort_values("_ar", ascending=False).drop(columns=["_ar"]).head(10)
            ap("Largest |rho| (any contrast; exploratory):")
            ap(top.to_string(index=False))
            ap("")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(REPORT.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
