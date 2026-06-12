"""
Extended integration: miRTarBase target union, MSigDB Hallmark ORA/GSEA,
Fisher overlaps (GSE188646 DE when present; HypoMap GEO tables with disclosure),
Enrichr aging-signature context, LIFU evidence layer,
multi-library ORA (global BH), bootstrap target-set stability, optional pseudobulk DE,
Stouffer combined p when two independent Fisher tests are available.
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import pandas as pd
from scipy.stats import combine_pvalues

from config import OUTPUT_DIR
import data_acquisition
import bootstrap_target_stability
import lifu_evidence_layer
import mirna_target_union
import msigdb_hallmark_gsea
import multi_library_ora
import overlap_fisher
import gse188646_sensitivity
import journal_tier_crossmodal
import journal_tier_crossmodal_celltype
import journal_tier_crosscohort_lability_replication
import journal_tier_niche_hypomap_external_validation
import journal_tier_niche_lability_localization
import journal_tier_niche_allen_aging_scrna_validation
import journal_tier_niche_merfish_spatial_validation
import journal_tier_niche_allen_aging_spatial_validation
import journal_tier_gse287202_v3_spot_audit
import journal_tier_allen_ish_marker_anatomy
import journal_tier_niche_small_strata_de
import journal_tier_crossmodal_meta_sensitivity
import journal_tier_negative_controls
import journal_tier_pathway_tf_coupling
import journal_tier_string_bridge
import replication_meta_two_cohorts
import sa_nsc_lifu_computational_suite
import tier2_analyses

HYPOMAP_DESEQ_URLS = [
    ("pomc_deseq", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208355/suppl/GSE208355_pomc_deseq.csv.gz"),
    ("agrp_deseq", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208355/suppl/GSE208355_agrp_deseq.csv.gz"),
    ("glp1r_deseq", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208355/suppl/GSE208355_glp1r_deseq.csv.gz"),
]


def _maybe_run_gse188646_stratified(log) -> None:
    col = os.environ.get("GSE188646_STRATUM_COL", "").strip()
    auto = os.environ.get("GSE188646_AUTO_STRATIFIED_PSEUDOBULK", "").strip().lower() in ("1", "true", "yes")
    man = OUTPUT_DIR / "gse188646_strata" / "manifest.csv"
    if not col and auto and not man.is_file():
        col = "seurat_clusters"
        log("Stratified pseudobulk: GSE188646_AUTO_STRATIFIED_PSEUDOBULK=1 → using seurat_clusters.")
    if not col:
        return
    rds_path = data_acquisition.resolve_gse188646_rds_path()
    if rds_path is None:
        log("Stratified pseudobulk: no GSE188646 RDS (set GSE188646_RDS or place canonical file under data/); skipped.")
        return
    r_script = Path(__file__).resolve().parent / "r" / "pseudobulk_stratified_edgeR_gse188646.R"
    if not r_script.is_file():
        log(f"Stratified R script not found: {r_script}")
        return
    cwd = Path(__file__).resolve().parent
    cmd = ["Rscript", str(r_script), str(rds_path), col]
    log(f"Running stratified pseudobulk (column {col}) — may take a long time...")
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=86400,
            encoding="utf-8",
            errors="replace",
        )
        if proc.returncode != 0:
            log(f"pseudobulk_stratified_edgeR_gse188646.R failed (exit {proc.returncode}).")
            if proc.stderr:
                log(proc.stderr[-4000:])
            return
        log("Stratified pseudobulk completed; see outputs/gse188646_strata/manifest.csv")
    except subprocess.TimeoutExpired:
        log("Stratified pseudobulk timed out after 86400s.")
    except FileNotFoundError:
        log("Rscript not found; stratified pseudobulk skipped.")


def _maybe_run_two_cohort_meta(log) -> None:
    """
    GSE87102 limma + DerSimonian–Laird meta vs GSE188646 DE.

    Default: ON when cohort1 DE exists. Opt-out: GSE188646_SKIP_COHORT2_META=1, or legacy
    GSE188646_RUN_COHORT2_META=0 (false/no).
    """
    if os.environ.get("GSE188646_SKIP_COHORT2_META", "").strip().lower() in ("1", "true", "yes"):
        log("Two-cohort meta: skipped (GSE188646_SKIP_COHORT2_META=1).")
        return
    if os.environ.get("GSE188646_RUN_COHORT2_META", "").strip().lower() in ("0", "false", "no"):
        log("Two-cohort meta: skipped (GSE188646_RUN_COHORT2_META=0).")
        return
    deg0 = OUTPUT_DIR / "gse188646_young_vs_aged_deg.csv"
    if not deg0.is_file():
        log("Two-cohort meta: outputs/gse188646_young_vs_aged_deg.csv missing; skipped.")
        return
    cwd = Path(__file__).resolve().parent
    r_c2 = cwd / "r" / "cohort2_gse87102_c57_hypothalamus_limma.R"
    if not r_c2.is_file():
        log(f"cohort2 R script missing: {r_c2}")
        return
    deg = OUTPUT_DIR / "gse188646_young_vs_aged_deg.csv"
    mtx = OUTPUT_DIR / "gse188646_pseudobulk_counts.mtx"
    need_se = True
    if deg.is_file():
        try:
            need_se = "se_logFC" not in pd.read_csv(deg, nrows=1).columns
        except Exception:
            need_se = True
    if need_se and mtx.is_file():
        r_reg = cwd / "r" / "regenerate_gse188646_deg_from_mtx.R"
        if r_reg.is_file():
            log(
                "Cohort1 DE missing se_logFC; regenerating from saved pseudobulk .mtx "
                "(r/regenerate_gse188646_deg_from_mtx.R) ..."
            )
            try:
                pr = subprocess.run(
                    ["Rscript", str(r_reg)],
                    check=False,
                    cwd=str(cwd),
                    capture_output=True,
                    text=True,
                    timeout=7200,
                    encoding="utf-8",
                    errors="replace",
                )
                if pr.returncode != 0:
                    log(f"regenerate_gse188646_deg_from_mtx.R failed (exit {pr.returncode}).")
                    if pr.stderr:
                        log(pr.stderr[-2500:])
                else:
                    log("Cohort1 DE updated with se_logFC.")
            except FileNotFoundError:
                log("Rscript not found; cannot regenerate cohort1 from mtx.")
            except subprocess.TimeoutExpired:
                log("regenerate_gse188646_deg_from_mtx.R timed out.")
        else:
            log(f"Missing {r_reg.name}; cannot backfill se_logFC from mtx.")
    elif need_se and not mtx.is_file():
        log("Cohort1 DE missing se_logFC and no pseudobulk .mtx; cannot run meta without RDS re-aggregation.")

    log("=== Cohort2 GSE87102 (C57 hypothalamus microarray) limma DE ===")
    try:
        proc = subprocess.run(
            ["Rscript", str(r_c2)],
            check=False,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=7200,
            encoding="utf-8",
            errors="replace",
        )
        if proc.returncode != 0:
            log(f"cohort2_gse87102_c57_hypothalamus_limma.R failed (exit {proc.returncode}).")
            if proc.stderr:
                log(proc.stderr[-4000:])
            return
        log("Cohort2 limma DE written to outputs/.")
    except FileNotFoundError:
        log("Rscript not found; cohort2/meta skipped.")
        return
    except subprocess.TimeoutExpired:
        log("cohort2 R timed out.")
        return

    log("=== Two-cohort random-effects meta (DL) on harmonised gene symbols ===")
    try:
        replication_meta_two_cohorts.run_two_cohort_meta(cwd, OUTPUT_DIR, log)
    except Exception as e:
        log(f"Two-cohort meta failed: {e}")


def _maybe_run_hypomap_reference_mapping(log) -> None:
    """Run HypoMap C185 expression-profile mapping when RDS + reference matrix exist."""
    p = data_acquisition.resolve_gse188646_rds_path()
    if p is None:
        return
    proj = Path(__file__).resolve().parent
    ref_default = proj / "data" / "references" / "hypomap_cellxgene_C185_named_mean_X_min200.csv"
    ref_env = os.environ.get("GSE188646_HYPOMAP_REF_EXPR_CSV", "").strip()
    ref = Path(ref_env) if ref_env else ref_default
    if not ref.is_file():
        log(
            "HypoMap reference mapping: no GSE188646_HYPOMAP_REF_EXPR_CSV and no default C185 matrix; skipped."
        )
        return
    out_csv = OUTPUT_DIR / "gse188646_hypomap_mapping" / "hypomap_custom_ref_spearman.csv"
    force = os.environ.get("GSE188646_FORCE_HYPOMAP_MAPPING", "").strip().lower() in ("1", "true", "yes")
    if out_csv.is_file() and not force:
        log("HypoMap custom-ref mapping exists; skipping R (set GSE188646_FORCE_HYPOMAP_MAPPING=1 to rerun).")
        return
    r_script = proj / "r" / "gse188646_hypomap_reference_mapping.R"
    if not r_script.is_file():
        log(f"HypoMap mapping R script not found: {r_script}")
        return
    env = os.environ.copy()
    env["GSE188646_HYPOMAP_REF_EXPR_CSV"] = str(ref.resolve())
    log(f"Running HypoMap reference mapping (ref={ref.name})...")
    try:
        proc = subprocess.run(
            ["Rscript", str(r_script), str(p)],
            check=False,
            cwd=str(proj),
            capture_output=True,
            text=True,
            timeout=3600,
            encoding="utf-8",
            errors="replace",
            env=env,
        )
        if proc.returncode != 0:
            log(f"HypoMap mapping R failed (exit {proc.returncode}): {proc.stderr[-2000:]}")
        else:
            log("HypoMap reference mapping completed.")
    except FileNotFoundError:
        log("Rscript not found; HypoMap mapping skipped.")
    except subprocess.TimeoutExpired:
        log("HypoMap mapping R timed out.")


def _maybe_run_gse188646_pseudobulk(log) -> None:
    """
    If GSE188646_RDS (or canonical data/GSE188646_*.RDS) points to the GEO Seurat object,
    run r/pseudobulk_edgeR_gse188646.R when DE output is missing or GSE188646_FORCE_PSEUDOBULK=1.
    Counts export: default skip (arg skip_counts). Set GSE188646_EXPORT_COUNTS=1 for .mtx/csv.
    """
    p = data_acquisition.resolve_gse188646_rds_path()
    if p is None:
        return
    deg = OUTPUT_DIR / "gse188646_young_vs_aged_deg.csv"
    force = os.environ.get("GSE188646_FORCE_PSEUDOBULK", "").strip().lower() in ("1", "true", "yes")
    if deg.exists() and not force:
        log(
            "GSE188646 young-vs-aged DE file exists; skipping R pseudobulk "
            "(set GSE188646_FORCE_PSEUDOBULK=1 to rerun)."
        )
        return
    r_script = Path(__file__).resolve().parent / "r" / "pseudobulk_edgeR_gse188646.R"
    if not r_script.is_file():
        log(f"R script not found: {r_script}")
        return
    cmd = ["Rscript", str(r_script), str(p)]
    export_counts = os.environ.get("GSE188646_EXPORT_COUNTS", "").strip().lower() in ("1", "true", "yes")
    if not export_counts:
        cmd.append("skip_counts")
    cwd = Path(__file__).resolve().parent
    log("Running pseudobulk edgeR on GSE188646 (this may take several minutes)...")
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=7200,
            encoding="utf-8",
            errors="replace",
        )
        if proc.returncode != 0:
            log(f"pseudobulk_edgeR_gse188646.R failed (exit {proc.returncode}).")
            if proc.stderr:
                log(proc.stderr[-4000:])
            return
        log("pseudobulk_edgeR_gse188646.R completed; outputs/gse188646_young_vs_aged_deg.csv updated.")
        if proc.stdout:
            tail = proc.stdout.strip().splitlines()[-5:]
            for line in tail:
                log(line)
    except subprocess.TimeoutExpired:
        log("pseudobulk Rscript timed out after 7200s.")
    except FileNotFoundError:
        log("Rscript not found on PATH; install R and retry.")


def _build_completeness_block(out_dir: Path, project_root: Path) -> list[str]:
    """Human-readable checklist for SA-oriented CSVs + figure bundle."""

    def row(label: str, p: Path, *, optional: bool = False) -> str:
        exists = p.is_file()
        if optional:
            if exists:
                st = "OK"
            else:
                st = "OPTIONAL absent"
        else:
            st = "OK" if exists else "MISSING"
        return f"  [{st}] {label}: {p.name}"

    lines = [
        "",
        "=== Science Advances bundle completeness ===",
        "  --- Prerequisite (run_pipeline.py) ---",
        row("MOESM htNSC vs astro miRNA summary", out_dir / "mirna_htnsc_astrocyte_summary.csv"),
        row("miRNA target union", out_dir / "mirna_target_union_genes.csv"),
        "  --- GSE188646 + cohort replication ---",
        row("GSE188646 pseudobulk DE", out_dir / "gse188646_young_vs_aged_deg.csv"),
        row("Cohort2 GSE87102 limma DE", out_dir / "cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv"),
        row("Two-cohort DL meta", out_dir / "exploratory_meta_DE_two_cohort_DL.csv"),
        "  --- Journal tier (primary computational spine) ---",
        row("Cross-modal per-gene burden", out_dir / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv"),
        row("Cross-modal summary JSON", out_dir / "exploratory_crossmodal_mirna_aging_summary.json"),
        row("Cross-modal gene-shuffle nulls", out_dir / "exploratory_crossmodal_permutation_rho_gene_shuffle.csv"),
        row("Cross-modal random-miRNA-set nulls", out_dir / "exploratory_crossmodal_random_mirna_set_rho_null.csv"),
        row("Stratified neg-control gene covariates", out_dir / "exploratory_negative_controls_gene_covariates.csv"),
        row("Stratified neg-control perm rhos", out_dir / "exploratory_negative_controls_stratified_perm_rhos.csv"),
        row("Stratified negative-controls summary", out_dir / "exploratory_negative_controls_summary.json"),
        row("Niche lability per-stratum table", out_dir / "exploratory_niche_lability_per_stratum.csv", optional=True),
        row("Niche lability localization summary", out_dir / "exploratory_niche_lability_localization_summary.json", optional=True),
        row("Pathway/TF coupling summary", out_dir / "exploratory_pathway_tf_coupling_summary.json", optional=True),
        row("STRING Piezo1 bridge summary", out_dir / "exploratory_string_piezo1_bridge_summary.json", optional=True),
        row("Cross-modal x meta sensitivity (genes)", out_dir / "exploratory_crossmodal_burden_vs_two_cohort_meta.csv"),
        row("Cross-modal x meta sensitivity (JSON)", out_dir / "exploratory_crossmodal_meta_cohort_sensitivity_summary.json"),
        row("Cross-cohort lability replication", out_dir / "exploratory_crosscohort_lability_replication_summary.json", optional=True),
        row("HypoMap external niche validation", out_dir / "exploratory_niche_hypomap_external_validation_summary.json", optional=True),
        "  --- SA niche / LIFU prior suite ---",
        row("SA curated GSEA", out_dir / "exploratory_sa_nsc_lifu_fgsea_curated_sets.csv"),
        row("SA suite summary JSON", out_dir / "exploratory_sa_nsc_lifu_suite_summary.json"),
        "  --- Tier-2 depth (optional; needs pseudobulk counts) ---",
        row("Pseudobulk counts .mtx", out_dir / "gse188646_pseudobulk_counts.mtx", optional=True),
        row("PROGENy pathway tests", out_dir / "progeny_pseudobulk_young_vs_aged_welch.csv", optional=True),
    ]
    fig_dir = out_dir / "figures" / "sa_bundle"
    man = project_root / "data" / "provenance" / "FIGURE_PANEL_MANIFEST.csv"
    expected = 5
    if man.is_file():
        try:
            import pandas as pd

            expected = len(pd.read_csv(man))
        except Exception:
            pass
    n_png = len(list(fig_dir.glob("*.png"))) if fig_dir.is_dir() else 0
    if n_png >= expected:
        pst = "OK"
    elif n_png:
        pst = "PARTIAL"
    else:
        pst = "MISSING"
    lines.append("  --- Figure bundle (manifest) ---")
    lines.append(f"  [{pst}] Figure PNGs in figures/sa_bundle/: {n_png} (manifest rows {expected})")
    lines.append(row("Figure manifest file", man))
    lines.append(
        row(
            "Pathway convergence Jaccard table",
            out_dir / "exploratory_pathway_convergence_jaccard.csv",
            optional=True,
        )
    )
    lines.append(
        row(
            "Pathway convergence figure PNG",
            out_dir / "figures" / "sa_bundle" / "fig_pathway_convergence.png",
            optional=True,
        )
    )
    lines.append(
        row(
            "Hallmark stress UpSet supplementary PNG",
            out_dir / "figures" / "sa_bundle" / "fig_supp_hallmark_stress_upset.png",
            optional=True,
        )
    )
    lines.append(
        row(
            "Strata cross-modal supplementary figures Word",
            out_dir / "manuscript" / "Supplementary_Figures_Strata_Crossmodal.docx",
            optional=True,
        )
    )
    lines.append(
        "  Opt-out env: GSE188646_SKIP_COHORT2_META=1 (no cohort2/meta); SKIP_SA_FIGURE_BUNDLE=1 (no PNG build)."
    )
    return lines


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "EXTENDED ANALYSIS — provenance",
        "- miRNA targets: Enrichr GMT miRTarBase_2017 (mmu-miR / mmu-let rows), experimentally reported interactions.",
        "- Hallmark ORA/GSEA: Enrichr GMT MSigDB_Hallmark_2020 (downloaded to data/cache).",
        "- Multi-library ORA: several mouse Enrichr libraries; global BH on pooled term P-values (exploratory).",
        "- Bootstrap: Jaccard stability of target union under resampling of top-ranked miRNAs (method sensitivity).",
        "- Aging context ORA: Enrichr library Aging_Perturbations_from_GEO_up (meta-signatures; not tissue-specific).",
        "- HypoMap: GEO GSE208355 supplementary *_deseq.csv.gz (IP/input contrasts; NOT hypothalamus young vs aged).",
        "- GSE188646: env GSE188646_RDS or canonical data/GSE188646_hypo.integrated.final.20210719.RDS runs "
        "r/pseudobulk_edgeR_gse188646.R (edgeR QLF; recommended). Optional: AUTO_FETCH_GSE188646_RDS=1, "
        "AUTO_FETCH_PUBLIC_DATA=1, AUTO_FETCH_HYPOMAP_H5AD=1 (see data/provenance/DATA_ACQUISITION_AND_QC_RATIONALE.txt).",
        "- Stouffer: combines Fisher P-values from GSE188646 and HypoMap overlaps only as a descriptive summary; "
        "contrasts are not exchangeable (different biology), so do not treat as a formal meta-analysis.",
        "- Two-cohort DE meta (default ON): GSE188646 pseudobulk + GSE87102 C57 hypothalamus microarray (limma), "
        "DerSimonian–Laird on logFC+SE — runs each extended pass when cohort1 DE exists; "
        "opt-out with GSE188646_SKIP_COHORT2_META=1 (see data/replication/META_TWO_COHORTS_README.txt).",
        "- Tier-2: R Bioconductor fgsea on the same Hallmark GMT (complements gseapy prerank); GenAge Mus musculus "
        "Fisher overlap (HAGR models_genes.zip; see data/provenance/GenAge_HAGR.txt); decoupler PROGENy and "
        "DoRothEA (mlm regulon activity) on pseudobulk when .mtx is present (export counts with GSE188646_EXPORT_COUNTS=1); "
        "pseudobulk Hallmark fgsea (multilevel when available), limma camera, limma fry; extra GEO sample tables.",
        "- LIFU: outputs/LIFU_evidence_layer.txt (literature layer; not statistically merged with MOESM).",
        "- Locked narrative + SAP: data/provenance/SCIENTIFIC_STORY_ONE_PAGE.txt, STATISTICAL_ANALYSIS_PLAN.txt, "
        "EXTERNAL_DATA_PROVENANCE_TABLE.csv; data/provenance/DATA_ACQUISITION_AND_QC_RATIONALE.txt; "
        "exploratory GSE188646 Fisher/camera/fgsea/fry/DoRothEA outputs use the exploratory_* prefix.",
        "- Journal-tier cross-modal (exploratory): miRNA target burden vs GSE188646 logFC with gene-shuffle "
        "and random-miRNA-set nulls (outputs/exploratory_crossmodal_*); burden vs two-cohort DL meta "
        "when exploratory_meta_DE_two_cohort_DL.csv exists (exploratory_crossmodal_burden_vs_two_cohort_meta.csv).",
        "- Journal-tier negative controls: se_logFC-strata + program-degree + GMT-wide targetability "
        "matched logFC shuffles (outputs/exploratory_negative_controls_*).",
        "- Cell-type / cluster cross-modal: outputs/gse188646_strata/* from R when "
        "GSE188646_STRATUM_COL or GSE188646_AUTO_STRATIFIED_PSEUDOBULK=1; Python summarizes "
        "exploratory_crossmodal_celltype_strata_*; niche lability localization "
        "(exploratory_niche_lability_*) tests whether |logFC| lability concentrates in "
        "third-ventricle tanycyte/NSC-niche clusters (marker-module labels); optional "
        "RUN_CROSSMODAL_STRATA_DIAGNOSTICS=1 runs tools/diagnostic_crossmodal_strata_figures.py; add "
        "RUN_CROSSMODAL_STRATA_SUPPLEMENTARY_DOCX=1 to also write Supplementary_Figures_Strata_Crossmodal.docx.",
        "- Pathway/TF coupling: exploratory_pathway_tf_coupling_summary.json and CSVs — burden vs "
        "PROGENy / DoRothEA Young–Aged deltas with gene-level burden permutation null.",
        "- STRING Piezo1 bridge: exploratory_string_piezo1_bridge_summary.json (STRING-db API; "
        "opt-out STRING_BRIDGE_OFFLINE=1).",
        "- SA NSC/LIFU suite (exploratory): curated niche/LIFU-prior gene sets vs aging DE (fgsea), "
        "pseudobulk module scores, Fisher vs miRNA targets, HypoMap niche flags (exploratory_sa_nsc_lifu_*).",
        "- Figure bundle (default ON): figures/build_sa_figure_bundle.py after validation; manifest "
        "data/provenance/FIGURE_PANEL_MANIFEST.csv; PNGs under outputs/figures/sa_bundle/; "
        "opt-out with SKIP_SA_FIGURE_BUNDLE=1.",
        "- Completeness checklist: outputs/SA_COMPLETENESS_CHECK.txt (mirrors final EXTENDED_REPORT section).",
        "",
    ]

    def log(msg: str):
        lines.append(msg)
        print(msg)

    log("\n=== Reference data: optional auto-fetch + GMT cache QC ===")
    try:
        data_acquisition.run_prefetch_hooks(log)
    except Exception as e:
        log(f"Prefetch hooks: {e}")
    try:
        data_acquisition.validate_existing_public_gmts(log)
    except Exception as e:
        log(f"GMT cache validation failed (delete corrupt files under data/cache if needed): {e}")
    data_acquisition.describe_rds_resolution(log)

    log("=== Optional: GSE188646 pseudobulk DE from RDS (env or canonical data/ path) ===")
    _maybe_run_gse188646_pseudobulk(log)
    _maybe_run_gse188646_stratified(log)
    _maybe_run_two_cohort_meta(log)

    log("=== Extended: miRNA -> miRTarBase (mmu) target union ===")
    mirna_csv = OUTPUT_DIR / "mirna_htnsc_astrocyte_summary.csv"
    if not mirna_csv.exists():
        raise FileNotFoundError("Run run_pipeline.py first to create mirna_htnsc_astrocyte_summary.csv")
    mirna_df_full = pd.read_csv(mirna_csv)
    top = mirna_target_union.load_top_htnsc_mirnas(mirna_csv, top_n=60, min_logfc=0.35)
    pd.Series(top, name="mirna").to_csv(OUTPUT_DIR / "top_htnsc_mirnas_for_targets.csv", index=False)
    gmt = mirna_target_union.ensure_mirtarbase_gmt()
    data_acquisition.validate_gmt_file(gmt, "miRTarBase_2017")
    hallmark_gmt_path = msigdb_hallmark_gsea.ensure_hallmark_gmt()
    data_acquisition.validate_gmt_file(hallmark_gmt_path, "MSigDB_Hallmark_2020")
    union, per_mir = mirna_target_union.union_targets_for_mirnas(top, gmt)
    mirna_target_union.save_target_maps(union, per_mir, OUTPUT_DIR)
    log(f"Union targets: {len(union)} genes from {len(per_mir)} miRNAs with miRTarBase entries")

    uni_mmu = mirna_target_union.mmu_gmt_gene_universe(gmt)
    ranked = mirna_target_union.build_ranked_genes_for_gsea(per_mir, uni_mmu)
    ranked.to_csv(OUTPUT_DIR / "mirna_target_multiplicity_rank.csv", index=False)

    log("\n=== MSigDB Hallmark ORA (mouse) on target union ===")
    ora_h = msigdb_hallmark_gsea.run_hallmark_ora(
        sorted(union),
        "miRTarBase_union_top60_htnsc_mirnas",
        background_genes=uni_mmu,
    )
    msigdb_hallmark_gsea.save_df(ora_h, OUTPUT_DIR / "enrichr_hallmark_ora_mirtarbase_union.csv")
    if "ora_backend" in ora_h.columns:
        ob = ora_h["ora_backend"].iloc[0]
        extra = ""
        if "ora_background_note" in ora_h.columns:
            n = str(ora_h["ora_background_note"].iloc[0])
            if n and n != "nan":
                extra = f" ({n})"
        log(f"Hallmark ORA backend: {ob}{extra}")
    log(ora_h.head(12).to_string(index=False))

    log("\n=== Multi-library Enrichr ORA + global BH (may take several minutes) ===")
    try:
        long_ora, sum_ora = multi_library_ora.run_multi_ora(
            sorted(union),
            libraries=None,
            organism="mouse",
            top_n_terms_per_lib=35,
        )
        multi_library_ora.save_multi_ora(long_ora, sum_ora, OUTPUT_DIR)
        log(sum_ora.to_string(index=False))
        if len(long_ora):
            topg = long_ora.sort_values("global_bh_q", na_position="last").head(15)
            log("\nTop terms by global BH q:\n" + topg.to_string(index=False))
    except Exception as e:
        log(f"Multi-library ORA skipped or partial failure: {e}")

    log("\n=== Bootstrap: target-union stability vs top-miRNA reference ===")
    try:
        bdf = bootstrap_target_stability.bootstrap_target_stability(
            mirna_df_full,
            mirna_col="mirna",
            score_col="logfc_htnsc_vs_astro",
            gmt_path=gmt,
            pool_size=80,
            n_draw=50,
            sample_size=35,
            seed=7,
        )
        bootstrap_target_stability.save_bootstrap(
            bdf, OUTPUT_DIR / "bootstrap_target_union_stability.csv"
        )
        summ = bdf.attrs.get("summary", {})
        log(str(summ))
    except Exception as e:
        log(f"Bootstrap stability skipped: {e}")

    log("\n=== Enrichr: Aging_Perturbations_from_GEO_up (context ORA) ===")
    try:
        ora_a = msigdb_hallmark_gsea.run_aging_geo_up_ora(
            sorted(union),
            "miRTarBase_union",
            multi_library_terms=OUTPUT_DIR / "enrichr_multi_library_ora_terms.csv",
        )
        msigdb_hallmark_gsea.save_df(ora_a, OUTPUT_DIR / "enrichr_aging_geo_up_ora_mirtarbase_union.csv")
        if len(ora_a):
            log(ora_a.head(15).to_string(index=False))
        else:
            log("Aging GEO ORA: empty results.")
    except Exception as e:
        log(f"Aging GEO ORA failed after retries / offline fallback: {e}")

    log("\n=== Preranked GSEA vs MSigDB Hallmark (target multiplicity) ===")
    try:
        tmp = ranked[["gene", "rank_metric"]].copy()
        gsea = msigdb_hallmark_gsea.run_prerank_hallmark_gsea(
            tmp,
            "gene",
            "rank_metric",
            OUTPUT_DIR / "hallmark_gsea_from_target_multiplicity",
        )
        msigdb_hallmark_gsea.save_df(gsea, OUTPUT_DIR / "gsea_hallmark_prerank_results.csv")
        if gsea is not None and len(gsea):
            log(gsea.head(12).to_string(index=False))
    except Exception as e:
        log(f"GSEA prerank skipped: {e}")

    fisher_p_gse = None
    fisher_p_hypo = None

    log("\n=== Fisher: targets vs GSE188646 DE (if R script output exists) ===")
    gse_deg = OUTPUT_DIR / "gse188646_young_vs_aged_deg.csv"
    if gse_deg.exists():
        ext = overlap_fisher._read_seurat_markers(gse_deg, padj_max=0.05)
        ext_in_u = {g.upper() for g in ext} & {g.upper() for g in uni_mmu}
        res = overlap_fisher.fisher_overlap(union, ext_in_u, uni_mmu)
        fisher_p_gse = res.get("fisher_two_sided_p")
        pd.Series(res).to_csv(OUTPUT_DIR / "fisher_targets_vs_gse188646_de.csv")
        log(str(res))
        if res.get("a_targets_and_external", 0) == 0 and res.get("c_external_only", 0) == 0:
            log(
                "Note: no GSE188646 DE genes passed padj<=0.05 in the overlap universe, so the Fisher "
                "test is degenerate (p=1). Consider relaxing padj for exploration only, or interpret as "
                "weak bulk-level separation after pseudobulk filtering."
            )
        gse188646_sensitivity.run_gse188646_fisher_sensitivity(union, uni_mmu, gse_deg, OUTPUT_DIR, log)
    else:
        log(
            "Skipped: outputs/gse188646_young_vs_aged_deg.csv not found. "
            "Place GSE188646_hypo.integrated.final.20210719.RDS under data/, set GSE188646_RDS, "
            "or AUTO_FETCH_GSE188646_RDS=1; then run r/pseudobulk_edgeR_gse188646.R / "
            "r/extract_gse188646_young_vs_aged.R (see data/provenance)."
        )

    log("\n=== HypoMap GSE208355 supplementary DESeq2 (NOT young-vs-aged; IP contrasts) ===")
    try:
        hypo_genes, hypo_meta = overlap_fisher.load_hypomap_de_union(
            HYPOMAP_DESEQ_URLS, padj_max=0.05, abs_lfc_min=0.5
        )
        hypo_meta.to_csv(OUTPUT_DIR / "hypomap_geo_deseq_filtered_concat.csv", index=False)
        uni_upper = {g.upper() for g in uni_mmu}
        hypo_upper = {g.upper() for g in hypo_genes}
        hypo_in_u = hypo_upper & uni_upper
        pd.Series(sorted(hypo_in_u), name="gene").to_csv(
            OUTPUT_DIR / "hypomap_de_genes_in_mirtarbase_universe.csv", index=False
        )
        res2 = overlap_fisher.fisher_overlap(union, hypo_in_u, uni_mmu)
        fisher_p_hypo = res2.get("fisher_two_sided_p")
        pd.Series(res2).to_csv(OUTPUT_DIR / "fisher_targets_vs_hypomap_geo_deseq_union.csv")
        log("HypoMap union genes: " + str(len(hypo_genes)))
        log(str(res2))
        log("See data/provenance/HypoMap_GSE208355.txt for interpretation limits.")
    except Exception as e:
        log(f"HypoMap download/overlap failed: {e}")

    log("\n=== Descriptive Stouffer combination of Fisher P-values (see provenance caveat) ===")
    ps = []
    labels = []
    if fisher_p_gse is not None:
        ps.append(float(fisher_p_gse))
        labels.append("gse188646_targets_vs_de")
    if fisher_p_hypo is not None:
        ps.append(float(fisher_p_hypo))
        labels.append("hypomap_targets_vs_ip_deseq")
    if len(ps) >= 2:
        degenerate = any(p >= 1.0 - 1e-12 or p <= 1e-12 for p in ps)
        if degenerate:
            meta_path = OUTPUT_DIR / "fisher_meta_stouffer_summary.csv"
            if meta_path.exists():
                meta_path.unlink(missing_ok=True)
            log(
                "Stouffer skipped: at least one Fisher p-value is degenerate (e.g. p=1 when the DE list "
                "is empty in the overlap universe). Combined p is not reported in that case."
            )
        else:
            comb = combine_pvalues(ps, method="stouffer")
            stat = float(comb.statistic)
            if not math.isfinite(stat):
                log("Stouffer produced a non-finite statistic; skipping CSV write.")
            else:
                meta_row = {
                    "n_tests": len(ps),
                    "labels": ";".join(labels),
                    "stouffer_statistic": stat,
                    "stouffer_combined_p": float(comb.pvalue),
                    "note": "HypoMap and GSE188646 are different contrasts; combined p is exploratory only.",
                }
                pd.DataFrame([meta_row]).to_csv(OUTPUT_DIR / "fisher_meta_stouffer_summary.csv", index=False)
                log(str(meta_row))
    else:
        log("Skipped: need both GSE188646 DE Fisher and HypoMap Fisher for Stouffer summary.")

    log("\n=== Journal-tier cross-modal coupling (exploratory) ===")
    try:
        journal_tier_crossmodal.run_crossmodal_mirna_aging(OUTPUT_DIR, log, n_perm=1000, n_mirna_draws=200)
    except Exception as e:
        log(f"Cross-modal coupling failed: {e}")
    try:
        journal_tier_crossmodal_meta_sensitivity.run_crossmodal_meta_sensitivity(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Cross-modal meta sensitivity failed: {e}")
    try:
        journal_tier_crosscohort_lability_replication.run_crosscohort_lability_replication(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Cross-cohort lability replication failed: {e}")
    try:
        journal_tier_negative_controls.run_expression_degree_negative_controls(OUTPUT_DIR, log, n_perm=800)
    except Exception as e:
        log(f"Negative-control pass failed: {e}")
    try:
        sa_nsc_lifu_computational_suite.run_suite(OUTPUT_DIR, log)
    except Exception as e:
        log(f"SA NSC/LIFU computational suite failed: {e}")

    proj = Path(__file__).resolve().parent
    tier2_analyses.run_tier2(
        log,
        project_root=proj,
        out_dir=OUTPUT_DIR,
        union=union,
        uni_mmu=uni_mmu,
        ranked_csv=OUTPUT_DIR / "mirna_target_multiplicity_rank.csv",
    )

    try:
        journal_tier_crossmodal_celltype.run_celltype_strata_crossmodal(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Cell-type / stratum cross-modal failed: {e}")
    try:
        journal_tier_niche_lability_localization.run_niche_lability_localization(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Niche lability localization failed: {e}")
    _maybe_run_hypomap_reference_mapping(log)
    try:
        journal_tier_niche_hypomap_external_validation.run_niche_hypomap_external_validation(OUTPUT_DIR, log)
    except Exception as e:
        log(f"HypoMap external niche validation failed: {e}")
    try:
        journal_tier_niche_allen_aging_scrna_validation.run_allen_aging_scrna_validation(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Allen Jin aging scRNA validation failed: {e}")
    try:
        journal_tier_niche_merfish_spatial_validation.run_merfish_spatial_validation(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Allen MERFISH spatial validation failed: {e}")
    try:
        journal_tier_niche_allen_aging_spatial_validation.run_allen_aging_spatial_validation(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Allen RSTE3 BIL spatial validation failed: {e}")
    try:
        journal_tier_gse287202_v3_spot_audit.run_gse287202_v3_spot_audit(OUTPUT_DIR, log)
    except Exception as e:
        log(f"GSE287202 V3 spot audit failed: {e}")
    try:
        journal_tier_allen_ish_marker_anatomy.run_allen_ish_marker_anatomy(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Allen ISH marker anatomy failed: {e}")
    try:
        journal_tier_niche_small_strata_de.run_small_niche_strata_de(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Small niche strata DE (48/50/57) failed: {e}")
    if os.environ.get("RUN_CROSSMODAL_STRATA_DIAGNOSTICS", "").strip().lower() in ("1", "true", "yes"):
        diag = proj / "tools" / "diagnostic_crossmodal_strata_figures.py"
        if diag.is_file() and (OUTPUT_DIR / "exploratory_crossmodal_celltype_strata_summary.csv").is_file():
            log("\n=== Cross-modal strata diagnostics (figures + annotated table) ===")
            try:
                cmd = [sys.executable, str(diag), "--outputs-dir", str(OUTPUT_DIR.resolve())]
                if os.environ.get("RUN_CROSSMODAL_STRATA_SUPPLEMENTARY_DOCX", "").strip().lower() in (
                    "1",
                    "true",
                    "yes",
                ):
                    cmd.append("--write-supplementary-docx")
                dr = subprocess.run(
                    cmd,
                    cwd=str(proj),
                    capture_output=True,
                    text=True,
                    timeout=300,
                    encoding="utf-8",
                    errors="replace",
                )
                if dr.stdout:
                    for ln in dr.stdout.strip().splitlines():
                        log(ln)
                if dr.returncode != 0:
                    log(f"diagnostic_crossmodal_strata_figures.py exit {dr.returncode}")
                    if dr.stderr:
                        log(dr.stderr[-2000:])
            except Exception as e:
                log(f"Strata diagnostics skipped: {e}")
    try:
        journal_tier_pathway_tf_coupling.run_pathway_tf_coupling(OUTPUT_DIR, log)
    except Exception as e:
        log(f"Pathway/TF coupling failed: {e}")
    try:
        journal_tier_string_bridge.run_string_piezo_bridge(OUTPUT_DIR, log)
    except Exception as e:
        log(f"STRING Piezo1 bridge failed: {e}")

    val_script = proj / "tools" / "validation_gse188646_priority_report.py"
    if val_script.is_file() and (OUTPUT_DIR / "gse188646_young_vs_aged_deg.csv").is_file():
        log("\n=== GSE188646 validation + prioritisation report ===")
        try:
            vr = subprocess.run(
                [sys.executable, str(val_script)],
                cwd=str(proj),
                capture_output=True,
                text=True,
                timeout=600,
                encoding="utf-8",
                errors="replace",
            )
            if vr.returncode == 0 and vr.stdout.strip():
                log(vr.stdout.strip().splitlines()[-1])
            elif vr.returncode != 0:
                log(f"validation_gse188646_priority_report.py exit {vr.returncode}")
                if vr.stderr:
                    log(vr.stderr[-2000:])
        except Exception as e:
            log(f"Validation report skipped: {e}")

    for tool_name in (
        "prepare_zenodo_deposit_manifest.py",
        "build_nature_portfolio_reporting_checklist.py",
        "build_zenodo_upload_bundle.py",
    ):
        tpath = proj / "tools" / tool_name
        if tpath.is_file():
            try:
                subprocess.run(
                    [sys.executable, str(tpath), "--outputs-dir", str(OUTPUT_DIR.resolve())],
                    cwd=str(proj),
                    check=False,
                    timeout=120,
                )
            except Exception as e:
                log(f"{tool_name} skipped: {e}")

    ms_script = proj / "tools" / "build_science_advances_manuscript_docx.py"
    if ms_script.is_file() and os.environ.get("SKIP_MANUSCRIPT_BUILD", "").strip().lower() not in (
        "1",
        "true",
        "yes",
    ):
        log("\n=== Manuscript DOCX build ===")
        try:
            mr = subprocess.run(
                [
                    sys.executable,
                    str(ms_script),
                    "--outputs-dir",
                    str(OUTPUT_DIR.resolve()),
                ],
                cwd=str(proj),
                capture_output=True,
                text=True,
                timeout=600,
                encoding="utf-8",
                errors="replace",
            )
            if mr.stdout:
                log(mr.stdout.strip().splitlines()[-1])
            if mr.returncode != 0 and mr.stderr:
                log(mr.stderr[-1500:])
        except Exception as e:
            log(f"Manuscript build skipped: {e}")

    if os.environ.get("SKIP_SA_FIGURE_BUNDLE", "").strip().lower() not in ("1", "true", "yes"):
        fig_script = proj / "figures" / "build_sa_figure_bundle.py"
        log("\n=== Science Advances figure bundle (manifest PNGs) ===")
        if fig_script.is_file():
            try:
                fr = subprocess.run(
                    [sys.executable, str(fig_script), "--outputs-dir", str(OUTPUT_DIR.resolve())],
                    cwd=str(proj),
                    capture_output=True,
                    text=True,
                    timeout=600,
                    encoding="utf-8",
                    errors="replace",
                )
                if fr.stdout:
                    for ln in fr.stdout.strip().splitlines()[-25:]:
                        log(ln)
                if fr.returncode != 0:
                    log(f"build_sa_figure_bundle.py exit {fr.returncode}")
                    if fr.stderr:
                        log(fr.stderr[-2000:])
            except Exception as e:
                log(f"Figure bundle skipped: {e}")
            pdf_script = proj / "tools" / "export_figure_bundle_pdf.py"
            if pdf_script.is_file():
                try:
                    pr = subprocess.run(
                        [sys.executable, str(pdf_script), "--outputs-dir", str(OUTPUT_DIR.resolve())],
                        cwd=str(proj),
                        capture_output=True,
                        text=True,
                        timeout=300,
                        encoding="utf-8",
                        errors="replace",
                    )
                    if pr.stdout:
                        for ln in pr.stdout.strip().splitlines()[-5:]:
                            log(ln)
                    if pr.returncode != 0 and pr.stderr:
                        log(pr.stderr[-1500:])
                except Exception as e:
                    log(f"Figure PDF export skipped: {e}")
        else:
            log(f"Missing {fig_script}")
    else:
        log("Figure bundle skipped (SKIP_SA_FIGURE_BUNDLE=1).")

    comp = _build_completeness_block(OUTPUT_DIR, proj)
    for line in comp:
        log(line)
    try:
        (OUTPUT_DIR / "SA_COMPLETENESS_CHECK.txt").write_text("\n".join(comp), encoding="utf-8")
    except Exception as e:
        log(f"Could not write SA_COMPLETENESS_CHECK.txt: {e}")

    p = lifu_evidence_layer.write_lifu_layer()
    log(f"\nWrote {p.name}")

    rep = OUTPUT_DIR / "EXTENDED_REPORT.txt"
    rep.write_text("\n".join(lines), encoding="utf-8")
    log(f"\nWrote {rep}")


if __name__ == "__main__":
    main()
