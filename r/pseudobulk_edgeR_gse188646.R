# Pseudobulk differential expression (recommended vs per-nucleus Wilcox for multi-animal designs)
# Requirements: Seurat, edgeR, dplyr, Matrix
#
# Usage (from feasibility_study/ directory):
#   Rscript r/pseudobulk_edgeR_gse188646.R /path/to/GSE188646_hypo.integrated.final.20210719.RDS
# Optional 2nd arg skip_counts|1|true|yes omits large count matrix export (DE still computed).
#
# Outputs:
#   outputs/gse188646_pseudobulk_counts.mtx (+ row/col name CSVs) or .csv if dense
#   outputs/gse188646_pseudobulk_metadata.csv
#   outputs/gse188646_young_vs_aged_deg.csv — edgeR QLF (gene, logFC, se_logFC, p_val, p_val_adj)
#
# Memory: set Sys.setenv(GSE188646_CELL_AGG_CHUNK = "60000") (or similar) before Rscript to sum
# counts per orig.ident in cell blocks (same result as one-shot AggregateExpression).
#
# Method: AggregateExpression sums counts per biological replicate (orig.ident),
#         then edgeR glmQLFTest on design ~ age (Aged vs Young).

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(dplyr)
  library(Matrix)
})

argv <- commandArgs(trailingOnly = FALSE)
fa <- grep("^--file=", argv, value = TRUE)
script_dir <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
source(file.path(script_dir, "gse188646_aggregate_counts.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript r/pseudobulk_edgeR_gse188646.R <seurat.rds>")
}
rds <- args[[1]]
if (!file.exists(rds)) stop("RDS not found: ", rds)

obj <- readRDS(rds)

guess_age <- function(md) {
  if ("age_bin" %in% colnames(md)) return(as.character(md$age_bin))
  if ("age_group" %in% colnames(md)) return(as.character(md$age_group))
  # GSE188646 integrated object: replicate labels in orig.ident (e.g. Young_1, Aged_3)
  if ("orig.ident" %in% colnames(md)) {
    oi <- as.character(md$orig.ident)
    ag <- rep(NA_character_, length(oi))
    ag[grepl("(?i)aged|old", oi, perl = TRUE)] <- "Aged"
    ag[grepl("(?i)young", oi, perl = TRUE)] <- "Young"
    if (!all(is.na(ag))) return(ag)
  }
  # GEO-style sample titles (when present on metadata)
  if ("title" %in% colnames(md)) {
    v <- as.character(md$title)
    ag <- rep(NA_character_, length(v))
    ag[grepl("(?i)aged|old", v, perl = TRUE)] <- "Aged"
    ag[grepl("(?i)young", v, perl = TRUE)] <- "Young"
    if (!all(is.na(ag))) return(ag)
  }
  cand <- intersect(c("orig.ident", "sample", "Sample"), colnames(md))
  if (length(cand) == 0) return(rep(NA_character_, nrow(md)))
  v <- as.character(md[[cand[[1]]]])
  ag <- rep(NA_character_, length(v))
  ag[grepl("(?i)young", v, perl = TRUE)] <- "Young"
  ag[grepl("(?i)aged|old", v, perl = TRUE)] <- "Aged"
  ag
}

md <- obj@meta.data
age <- guess_age(md)
if (all(is.na(age))) {
  stop("Could not infer Young/Aged. Add column age_bin to metadata or fix guess_age().")
}
obj$age_bin <- age
obj <- subset(obj, subset = !is.na(age_bin))

if ("RNA" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "RNA"
} else {
  warning("RNA assay missing; using ", DefaultAssay(obj))
}

# Pseudobulk = sum counts per orig.ident (GEO: one pooled hypothalamus pair per GSM title)
if (!"orig.ident" %in% colnames(obj@meta.data)) {
  stop("orig.ident required for pseudobulk replicate labels.")
}

# Do not set Idents(orig.ident): Seurat replaces '_' with '-' in identity names,
# which then mismatches metadata orig.ident (e.g. Aged_1 vs Aged-1) and breaks
# downstream alignment with pseudobulk matrix column names.
chunk_raw <- Sys.getenv("GSE188646_CELL_AGG_CHUNK", unset = "0")
chunk_size <- suppressWarnings(as.integer(chunk_raw))
if (length(chunk_size) != 1L || is.na(chunk_size)) {
  chunk_size <- 0L
}
counts <- aggregate_pseudobulk_counts(obj, chunk_size)
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "CsparseMatrix")
}

# Map each orig.ident to age
meta <- unique(obj@meta.data[, c("orig.ident", "age_bin"), drop = FALSE])
rownames(meta) <- meta$orig.ident
common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
meta$age_bin <- factor(meta$age_bin, levels = c("Young", "Aged"))

out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

skip_counts <- length(args) >= 2 && tolower(args[[2]]) %in% c("skip_counts", "1", "true", "yes")
if (skip_counts) {
  message("Skipping large pseudobulk counts export (arg 2 = skip_counts).")
} else if (inherits(counts, "sparseMatrix")) {
  mtx_path <- file.path(out_dir, "gse188646_pseudobulk_counts.mtx")
  Matrix::writeMM(counts, mtx_path)
  utils::write.csv(
    data.frame(gene = rownames(counts)),
    file.path(out_dir, "gse188646_pseudobulk_counts_rownames.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(sample = colnames(counts)),
    file.path(out_dir, "gse188646_pseudobulk_counts_colnames.csv"),
    row.names = FALSE
  )
  message("Wrote sparse counts: ", mtx_path)
} else {
  utils::write.csv(as.matrix(counts), file.path(out_dir, "gse188646_pseudobulk_counts.csv"))
}
utils::write.csv(meta, file.path(out_dir, "gse188646_pseudobulk_metadata.csv"), row.names = FALSE)

y <- DGEList(counts = counts, group = meta$age_bin)
keep <- filterByExpr(y, group = meta$age_bin)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ age_bin, data = meta)
colnames(design) <- sub("^age_bin", "", colnames(design))

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
# coef 2 = Aged vs Young (Young is intercept)
qlf <- glmQLFTest(fit, coef = "Aged")
tt <- topTags(qlf, n = Inf, sort.by = "none")
tab <- as.data.frame(tt$table)
# se from quasi-F (df1=1): |logFC| / sqrt(F); for meta-analysis with cohort 2 (REPLICATION_AND_META_POLICY).
fcol <- if ("F" %in% colnames(tab)) "F" else if ("LR" %in% colnames(tab)) "LR" else NULL
if (is.null(fcol)) {
  stop("edgeR qlf table missing F/LR column; cannot derive se_logFC.")
}
se_lfc <- abs(tab$logFC) / sqrt(pmax(tab[[fcol]], .Machine$double.eps))
res <- data.frame(
  gene = rownames(tab),
  logFC = tab$logFC,
  se_logFC = as.numeric(se_lfc),
  p_val = tab$PValue,
  p_val_adj = tab$FDR,
  stringsAsFactors = FALSE
)

utils::write.csv(res, file.path(out_dir, "gse188646_young_vs_aged_deg.csv"), row.names = FALSE)
message("Wrote pseudobulk DE: ", normalizePath(file.path(out_dir, "gse188646_young_vs_aged_deg.csv")))
