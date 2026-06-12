# Extract young vs aged pseudobulk or nucleus-level DE from GSE188646 Seurat RDS.
# Prerequisites (install in R): Seurat (v4+), dplyr
#
# Usage:
#   Rscript r/extract_gse188646_young_vs_aged.R path/to/GSE188646_hypo.integrated.final.20210719.RDS
#
# Output (written relative to feasibility_study/):
#   outputs/gse188646_young_vs_aged_deg.csv
#
# Notes:
# - DefaultAssay is set to "RNA" for count-based DE when available.
# - FindMarkers is a starting point; for publication-grade DE across animals,
#   prefer pseudobulk edgeR/DESeq2 on age × replicate (muscat / milo / aggregateExpression).
# - If metadata layout differs from expectations, inspect colnames(obj@meta.data)
#   and edit the age-guessing block below.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[[1]] else {
  stop("Provide RDS path: Rscript r/extract_gse188646_young_vs_aged.R <path_to_seurat.rds>")
}

obj <- readRDS(rds)
md <- obj@meta.data

guess_age <- function(md) {
  if ("age_bin" %in% colnames(md)) return(as.character(md$age_bin))
  if ("age_group" %in% colnames(md)) return(as.character(md$age_group))
  if ("AgeGroup" %in% colnames(md)) return(as.character(md$AgeGroup))
  if ("condition" %in% colnames(md)) return(as.character(md$condition))
  # GSE188646 integrated object: replicate labels in orig.ident (e.g. Young_1, Aged_3)
  if ("orig.ident" %in% colnames(md)) {
    oi <- as.character(md$orig.ident)
    ag <- rep(NA_character_, length(oi))
    ag[grepl("(?i)aged|old", oi, perl = TRUE)] <- "Aged"
    ag[grepl("(?i)young", oi, perl = TRUE)] <- "Young"
    if (!all(is.na(ag))) return(ag)
  }
  if ("title" %in% colnames(md)) {
    v <- as.character(md$title)
    ag <- rep(NA_character_, length(v))
    ag[grepl("(?i)aged|old", v, perl = TRUE)] <- "Aged"
    ag[grepl("(?i)young", v, perl = TRUE)] <- "Young"
    if (!all(is.na(ag))) return(ag)
  }
  cand_cols <- intersect(
    c("sample", "Sample", "library", "Library", "sample_id"),
    colnames(md)
  )
  if (length(cand_cols) == 0) return(rep(NA_character_, nrow(md)))
  v <- as.character(md[[cand_cols[[1]]]])
  ag <- rep(NA_character_, length(v))
  ag[grepl("(?i)young", v, perl = TRUE)] <- "Young"
  ag[grepl("(?i)aged|old", v, perl = TRUE)] <- "Aged"
  ag
}

age_bin <- guess_age(md)
if (all(is.na(age_bin))) {
  stop("Could not infer Young/Aged labels from metadata. Inspect colnames(obj@meta.data).")
}

obj$age_bin <- age_bin
obj <- subset(obj, subset = !is.na(age_bin))

if ("RNA" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "RNA"
} else {
  warning("RNA assay not found; using default assay: ", DefaultAssay(obj))
}

Idents(obj) <- "age_bin"
markers <- FindMarkers(
  obj,
  ident.1 = "Aged",
  ident.2 = "Young",
  test.use = "wilcox",
  logfc.threshold = 0.1,
  min.pct = 0.05,
  verbose = FALSE
)
markers$gene <- rownames(markers)

out_dir <- file.path("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_csv <- file.path(out_dir, "gse188646_young_vs_aged_deg.csv")
write.csv(markers, out_csv, row.names = FALSE)
message("Wrote ", normalizePath(out_csv))
