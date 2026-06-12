Journal-tier computational add-ons (exploratory)
===============================================

Mirror of outputs/JOURNAL_TIER_COMPUTATIONAL_README.txt (outputs/ may be gitignored).

This folder’s automated analyses are documented in EXTENDED_REPORT.txt and
data/provenance/SCIENTIFIC_STORY_ONE_PAGE.txt. The following are additional
computational layers aimed at top-tier *hypothesis framing* — not substitutes
for independent experiments or contrast-matched replication.

1) Cross-modal miRNA target burden vs GSE188646 aging DE
   - Module: src/journal_tier_crossmodal.py (invoked from run_extended.py)
   - Outputs (exploratory_*):
       exploratory_crossmodal_gene_burden_vs_aging_logfc.csv
       exploratory_crossmodal_mirna_aging_summary.json
       exploratory_crossmodal_permutation_rho_gene_shuffle.csv
       exploratory_crossmodal_random_mirna_set_rho_null.csv
   - Cross-modal × two-cohort DL meta (sensitivity; requires exploratory_meta_DE_two_cohort_DL.csv):
       exploratory_crossmodal_burden_vs_two_cohort_meta.csv
       exploratory_crossmodal_meta_cohort_sensitivity_summary.json
       (Module: src/journal_tier_crossmodal_meta_sensitivity.py; after cross-modal in run_extended.py.)

2) Negative controls: precision / program degree / GMT-wide targetability
   - Module: src/journal_tier_negative_controls.py (invoked after cross-modal in run_extended.py)
   - Outputs:
       exploratory_negative_controls_gene_covariates.csv
         (per gene: se_logFC, program n_mirnas, program weighted_burden, GMT-wide MOESM-weighted indegree)
       exploratory_negative_controls_stratified_perm_rhos.csv
       exploratory_negative_controls_summary.json
         (includes completed_utc, methodology_note, permutation p-values)
   - Nulls: Spearman(weighted_burden, logFC) after shuffling logFC **within** strata:
       (a) se_logFC deciles only (precision / detection proxy),
       (b) capped program target count x se deciles,
       (c) GMT-wide weighted in-degree decile x se deciles (network exposure beyond the top program).

   Interpretation: MOESM htNSC program vs hypothalamus aging DE remains non-exchangeable; these
   nulls stress-test whether any burden–logFC alignment could be driven by DE precision, local
   program degree, or global miRTarBase targetability.

   Why not a full degree-preserving bipartite graph null?
   A strict null would redraw edges (or swap endpoints) while fixing each miRNA's out-degree
   and each gene's in-degree (configuration / stub-matching). That is appropriate when edges are
   exchangeable and the graph is small enough for many independent draws. Here, edges are
   literature-supported priors (not exchangeable), the bipartite structure is large, and we hold
   the observed program burden vector w fixed to test Spearman(w, logFC). The implemented approach
   is the standard omics compromise: discretize continuous summaries (se_logFC deciles; capped
   program in-degree; GMT-wide weighted in-degree deciles) and shuffle logFC within strata so
   marginal confounding by precision and "targetability" is absorbed without claiming a full
   generative graph null. Random miRNA-set resampling (cross-modal module) addresses program
   specificity; GMT-wide indegree strata address global network exposure.

2a) Cell-type / cluster–resolved burden vs aging pseudobulk DE
   - R: r/pseudobulk_stratified_edgeR_gse188646.R → outputs/gse188646_strata/manifest.csv and
     stratum_*_young_vs_aged_deg.csv (requires GSE188646_STRATUM_COL, e.g. seurat_clusters, or set
     GSE188646_AUTO_STRATIFIED_PSEUDOBULK=1 to default seurat_clusters when manifest is absent).
   - Python: src/journal_tier_crossmodal_celltype.py (after Tier-2 in run_extended.py)
   - Outputs: exploratory_crossmodal_celltype_strata_summary.csv (+ .json)
   - Null: gene-label shuffle of logFC within each stratum’s merged gene list (same spirit as global
     cross-modal, but restricted to genes tested in that stratum’s DE table).
   - Optional diagnostics: python tools/diagnostic_crossmodal_strata_figures.py (forest ρ plot,
     burden–logFC scatter grid, glial-only forest, BH across strata on annotated CSV), or set RUN_CROSSMODAL_STRATA_DIAGNOSTICS=1
     when running python run_extended.py to invoke the same tool automatically; RUN_CROSSMODAL_STRATA_SUPPLEMENTARY_DOCX=1
     also writes outputs/manuscript/Supplementary_Figures_Strata_Crossmodal.docx (requires python-docx).

2b) Pathway- and TF-activity–level coupling (burden vs decoupler deltas)
   - Python: src/journal_tier_pathway_tf_coupling.py (after Tier-2; needs pseudobulk counts export
     for PROGENy / DoRothEA tables).
   - Outputs: exploratory_pathway_tf_coupling_summary.json,
     exploratory_pathway_burden_vs_progeny_delta.csv,
     exploratory_tf_burden_vs_dorothea_delta.csv
   - Null: permute per-gene weighted burden across genes in the cross-modal burden table while
     holding Young–Aged activity deltas fixed, then recompute pathway/TF mean burdens.

2c) STRING graph bridge (Piezo1 / mechanistic priors vs high-burden union genes)
   - Python: src/journal_tier_string_bridge.py (STRING-db.org API; set STRING_BRIDGE_OFFLINE=1 to skip).
   - Output: exploratory_string_piezo1_bridge_summary.json
   - Metrics: edges between fixed mechanism seeds and high-burden union subset in the induced
     network; stress-gene median graph distance to Piezo1∪1-hop with random gene-set null.
   - Null (edges): uniform random gene sets without replacement from the network node pool when
     |pool| ≥ |union subset|; otherwise degree-bin–matched draws. If |union| exceeds |pool|, the
     evaluated union subset is trimmed to the top-burden genes that fit the pool (documented in JSON).

3) Standard extended stack (unchanged contract)
   - miRTarBase union, ORA/GSEA, Fisher vs GSE188646 / HypoMap, bootstrap
     stability, tier-2 fgsea / PROGENy / DoRothEA / GenAge, optional two-cohort
     meta — see EXTENDED_REPORT.txt after each run.

4) Science Advances framing suite — curated niche / LIFU-biophysics priors (exploratory)
   - Module: src/sa_nsc_lifu_computational_suite.py (end of run_extended.py)
   - GMT written to outputs/sa_nsc_lifu_curated_sets.gmt (intersection with DE gene universe)
   - Outputs (exploratory_sa_nsc_lifu_*):
       exploratory_sa_nsc_lifu_fgsea_curated_sets.csv
       exploratory_sa_nsc_lifu_fisher_targets_vs_curated_sets.csv
       exploratory_sa_nsc_lifu_pseudobulk_module_young_vs_aged.csv
       exploratory_sa_nsc_lifu_burden_logfc_by_niche_subset.csv
       exploratory_sa_nsc_lifu_hypomap_ref_flags.csv (optional; requires hypomap_custom_ref_spearman.csv)
       exploratory_sa_nsc_lifu_hypomap_rho_summary_by_niche_flag.csv (optional)
       exploratory_sa_nsc_lifu_cluster_niche_flags.csv (optional; requires cluster_putative_labels.csv)
       exploratory_sa_nsc_lifu_suite_summary.json (completed_utc, methodology_note)
   - Interpretation: small literature gene sets stress whether aging DE / pseudobulk / cross-modal
     burden aligns with named hypothalamic niche and mechanosensitive / inflammatory priors;
     not merged with LIFU evidence layer statistically. Full figure mapping:
     data/provenance/FRAMING_NSC_HYPOTHALAMUS_LIFU_AGING_SCIENCE_ADVANCES.txt.

5) Figure bundle (manifest-driven PNGs; default at end of run_extended.py)
   - Manifest: data/provenance/FIGURE_PANEL_MANIFEST.csv
   - Builder: figures/build_sa_figure_bundle.py → outputs/figures/sa_bundle/*.png
   - Opt-out: SKIP_SA_FIGURE_BUNDLE=1
   - Machine checklist: outputs/SA_COMPLETENESS_CHECK.txt (mirrors EXTENDED_REPORT tail) asserts
     journal-tier files above (including random-miRNA-set null CSV and stratified neg-control
     covariates + perm-rho CSVs), cohort2 meta, SA suite, optional Tier-2, and manifest PNG count.

Re-run full stack from feasibility_study/:
  python run_extended.py

Optional env (Windows PowerShell examples):
  $env:GSE188646_RDS = "data\\GSE188646_hypo.integrated.final.20210719.RDS"
  $env:GSE188646_FORCE_PSEUDOBULK = "0"
  $env:GSE188646_SKIP_COHORT2_META = "1"
  $env:SKIP_SA_FIGURE_BUNDLE = "1"
  $env:GSE188646_STRATUM_COL = "seurat_clusters"
  # or auto-default seurat_clusters when manifest missing:
  $env:GSE188646_AUTO_STRATIFIED_PSEUDOBULK = "1"
  $env:STRING_BRIDGE_OFFLINE = "1"
