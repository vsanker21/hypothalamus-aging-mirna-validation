Two-cohort differential-expression meta-analysis (pre-specified in REPLICATION_AND_META_POLICY.txt)
================================================================================================

Cohort 1 (outputs/gse188646_young_vs_aged_deg.csv)
  GSE188646 integrated snRNA-seq object, pseudobulk edgeR QLF (Aged vs Young), 4 vs 4 female mice
  (GEO overall design). Columns include se_logFC derived from the quasi-likelihood F statistic.

Cohort 2 (outputs/cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv)
  GEO GSE87102: Agilent whole-mouse microarray, C57BL/6 whole hypothalamus, young (2–3 mo) vs old
  (22–24 mo) male replicates (sample titles on GEO). limma moderated t; se_logFC = |logFC|/|t|.

Contrast alignment
  Both files use positive logFC = higher expression in aged / old relative to young.

Honest design limits (do not over-claim)
  Sex differs (GSE188646 female; GSE87102 male in the selected C57 hypothalamus subset).
  Assay differs (snRNA pseudobulk vs bulk microarray). tau^2 from DerSimonian–Laird with k=2
  studies is unstable — use direction concordance, rank of |beta_DL|, and FDR only as a screen.

How to run (PowerShell examples)
  Two-cohort meta runs by default at the start of python run_extended.py when
  outputs/gse188646_young_vs_aged_deg.csv exists (requires R + GEOquery + limma for cohort2).

  Opt-out (faster CI or no R):
    $env:GSE188646_SKIP_COHORT2_META = "1"
    python run_extended.py

  Legacy explicit opt-out (same effect):
    $env:GSE188646_RUN_COHORT2_META = "0"
    python run_extended.py

  To refresh cohort1 DE with se_logFC before meta:
    $env:GSE188646_FORCE_PSEUDOBULK = "1"   # refresh cohort1 DE with se_logFC if CSV already exists
  $env:GSE188646_CELL_AGG_CHUNK = "60000" # optional: chunk AggregateExpression by cells (memory)
  $env:GSE188646_STRATUM_COL = "seurat_clusters"  # optional: per-cluster pseudobulk DE (many strata)
  python run_extended.py

Outputs
  outputs/exploratory_meta_DE_two_cohort_DL.csv
  outputs/exploratory_META_TWO_COHORT_DL_README.txt

Cross-modal × meta sensitivity (same miRNA burden vector vs pooled aging beta; exploratory)
  Requires exploratory_crossmodal_gene_burden_vs_aging_logfc.csv from run_extended journal-tier pass
  plus the meta CSV above. Writes:
  outputs/exploratory_crossmodal_burden_vs_two_cohort_meta.csv
  outputs/exploratory_crossmodal_meta_cohort_sensitivity_summary.json
  (Invoked automatically in run_extended.py after the cross-modal block whenever the meta CSV exists.)

GEO series matrix for GSE87102 is cached under data/cache/ after first GEOquery download.

If GSE188646_RDS is unavailable but outputs/gse188646_pseudobulk_counts.mtx exists, run:
  Rscript r/regenerate_gse188646_deg_from_mtx.R
to refresh cohort1 DE with se_logFC (same edgeR QLF as pseudobulk_edgeR after aggregation).
run_extended.py runs cohort2 + meta by default when cohort1 DE exists (unless GSE188646_SKIP_COHORT2_META=1).
