# Stratified pseudobulk edgeR: one DE table per metadata stratum (e.g. seurat_clusters).
# Use for cell-type / cluster–resolved young vs aged within the same GSE188646 animals.
#
# Usage (from feasibility_study/):
#   Rscript r/pseudobulk_stratified_edgeR_gse188646.R <seurat.rds> <stratum_colname>
#
# Env:
#   GSE188646_CELL_AGG_CHUNK — same as pseudobulk_edgeR_gse188646.R (optional)
#   GSE188646_STRATUM_MIN_CELLS — default 400
#   GSE188646_STRATUM_MIN_REP — min distinct orig.ident per age_bin (default 2)
#
# Outputs under outputs/gse188646_strata/:
#   manifest.csv — stratum, n_cells, n_young_rep, n_aged_rep, deg_path
#   stratum_<id>_young_vs_aged_deg.csv — per-stratum edgeR (same columns as main pseudobulk)

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
if (length(args) < 2) {
  stop("Usage: Rscript r/pseudobulk_stratified_edgeR_gse188646.R <seurat.rds> <stratum_col>")
}
rds <- args[[1]]
stratum_col <- args[[2]]
if (!file.exists(rds)) stop("RDS not found: ", rds)

chunk_raw <- Sys.getenv("GSE188646_CELL_AGG_CHUNK", unset = "0")
chunk_size <- suppressWarnings(as.integer(chunk_raw))
if (length(chunk_size) != 1L || is.na(chunk_size)) chunk_size <- 0L

min_cells <- suppressWarnings(as.integer(Sys.getenv("GSE188646_STRATUM_MIN_CELLS", unset = "400")))
if (is.na(min_cells)) min_cells <- 400L
min_genes <- suppressWarnings(as.integer(Sys.getenv("GSE188646_STRATUM_MIN_GENES", unset = "100")))
if (is.na(min_genes)) min_genes <- 100L
min_rep <- suppressWarnings(as.integer(Sys.getenv("GSE188646_STRATUM_MIN_REP", unset = "2")))
if (is.na(min_rep)) min_rep <- 2L
only_raw <- Sys.getenv("GSE188646_STRATUM_ONLY", unset = "")
only_strata <- character(0)
if (nzchar(only_raw)) {
  only_strata <- trimws(strsplit(only_raw, ",", fixed = TRUE)[[1L]])
  only_strata <- only_strata[nzchar(only_strata)]
}
append_manifest <- identical(tolower(Sys.getenv("GSE188646_STRATUM_APPEND_MANIFEST", unset = "0")), "true") ||
  Sys.getenv("GSE188646_STRATUM_APPEND_MANIFEST", unset = "0") %in% c("1", "yes")

guess_age <- function(md) {
  if ("age_bin" %in% colnames(md)) return(as.character(md$age_bin))
  if ("age_group" %in% colnames(md)) return(as.character(md$age_group))
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
  cand <- intersect(c("orig.ident", "sample", "Sample"), colnames(md))
  if (length(cand) == 0) return(rep(NA_character_, nrow(md)))
  v <- as.character(md[[cand[[1]]]])
  ag <- rep(NA_character_, length(v))
  ag[grepl("(?i)young", v, perl = TRUE)] <- "Young"
  ag[grepl("(?i)aged|old", v, perl = TRUE)] <- "Aged"
  ag
}

obj <- readRDS(rds)
md <- obj@meta.data
if (!stratum_col %in% colnames(md)) {
  stop(
    "Column not found in metadata: ", stratum_col,
    "\nAvailable columns (first 50): ",
    paste(head(colnames(md), 50L), collapse = ", "),
    "\nRe-run with a valid stratum column (e.g. seurat_clusters, RNA_snn_res.0.5).",
    call. = FALSE
  )
}
if (!"orig.ident" %in% colnames(md)) {
  stop("orig.ident required.")
}

age <- guess_age(md)
if (all(is.na(age))) stop("Could not infer Young/Aged.")
obj$age_bin <- age
obj <- subset(obj, subset = !is.na(age_bin))

if ("RNA" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "RNA"
}

strata <- sort(unique(as.character(obj@meta.data[[stratum_col]])))
strata <- strata[!is.na(strata) & nzchar(as.character(strata))]
if (length(only_strata) > 0L) {
  strata <- intersect(strata, only_strata)
  if (length(strata) == 0L) {
    stop("GSE188646_STRATUM_ONLY requested strata not found in ", stratum_col, ": ", only_raw)
  }
}

out_base <- file.path("outputs", "gse188646_strata")
if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)

manifest_rows <- list()

sanitize <- function(s) {
  gsub("[^A-Za-z0-9._-]+", "_", s, perl = TRUE)
}

for (st in strata) {
  cells_keep <- colnames(obj)[as.character(obj@meta.data[[stratum_col]]) == st]
  if (length(cells_keep) < min_cells) next
  sub <- subset(obj, cells = cells_keep)
  if (ncol(sub) < min_cells) next
  umeta <- unique(sub@meta.data[, c("orig.ident", "age_bin"), drop = FALSE])
  n_y <- sum(umeta$age_bin == "Young")
  n_a <- sum(umeta$age_bin == "Aged")
  if (n_y < min_rep || n_a < min_rep) next

  message("Stratum ", stratum_col, "=", st, " (", ncol(sub), " cells, ", nrow(umeta), " animals) ...")
  counts <- aggregate_pseudobulk_counts(sub, chunk_size)
  if (!inherits(counts, "sparseMatrix")) {
    counts <- as(counts, "CsparseMatrix")
  }
  meta <- unique(sub@meta.data[, c("orig.ident", "age_bin"), drop = FALSE])
  rownames(meta) <- meta$orig.ident
  common <- intersect(colnames(counts), rownames(meta))
  if (length(common) < min_rep * 2L) next
  counts <- counts[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  meta$age_bin <- factor(meta$age_bin, levels = c("Young", "Aged"))

  y <- DGEList(counts = counts, group = meta$age_bin)
  keep <- filterByExpr(y, group = meta$age_bin)
  if (sum(keep) < min_genes) next
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ age_bin, data = meta)
  colnames(design) <- sub("^age_bin", "", colnames(design))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "Aged")
  tt <- topTags(qlf, n = Inf, sort.by = "none")
  tab <- as.data.frame(tt$table)
  fcol <- if ("F" %in% colnames(tab)) "F" else if ("LR" %in% colnames(tab)) "LR" else NULL
  if (is.null(fcol)) next
  se_lfc <- abs(tab$logFC) / sqrt(pmax(tab[[fcol]], .Machine$double.eps))
  res <- data.frame(
    gene = rownames(tab),
    logFC = tab$logFC,
    se_logFC = as.numeric(se_lfc),
    p_val = tab$PValue,
    p_val_adj = tab$FDR,
    stringsAsFactors = FALSE
  )

  fn <- file.path(out_base, paste0("stratum_", sanitize(st), "_young_vs_aged_deg.csv"))
  utils::write.csv(res, fn, row.names = FALSE)
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    stratum_col = stratum_col,
    stratum = st,
    n_cells = ncol(sub),
    n_young_rep = n_y,
    n_aged_rep = n_a,
    deg_csv = basename(fn),
    stringsAsFactors = FALSE
  )
}

if (length(manifest_rows) == 0L) {
  stop("No strata passed filters (cells / replicates). Relax GSE188646_STRATUM_MIN_* or check column.")
}
man <- do.call(rbind, manifest_rows)
man_path <- file.path(out_base, "manifest.csv")
if (append_manifest && file.exists(man_path)) {
  old <- utils::read.csv(man_path, stringsAsFactors = FALSE)
  man <- rbind(old, man)
  man <- man[!duplicated(man[, c("stratum_col", "stratum")]), , drop = FALSE]
}
utils::write.csv(man, man_path, row.names = FALSE)
message("Wrote ", nrow(man), " stratified DE tables under ", normalizePath(out_base, winslash = "/", mustWork = FALSE))
