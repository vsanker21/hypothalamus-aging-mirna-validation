# Regenerate outputs/gse188646_young_vs_aged_deg.csv from saved pseudobulk counts + metadata
# (same edgeR QLF path as pseudobulk_edgeR_gse188646.R after aggregation). Use when RDS is
# unavailable but outputs/gse188646_pseudobulk_counts.mtx and metadata CSV already exist.
#
# Usage (from feasibility_study/):
#   Rscript r/regenerate_gse188646_deg_from_mtx.R
# Or custom paths:
#   Rscript r/regenerate_gse188646_deg_from_mtx.R <mtx> <rownames.csv> <colnames.csv> <metadata.csv> <out_deg.csv>

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- "outputs"
if (length(args) >= 5) {
  mtx_path <- args[[1]]
  rn_path <- args[[2]]
  cn_path <- args[[3]]
  meta_path <- args[[4]]
  out_deg <- args[[5]]
} else {
  mtx_path <- file.path(out_dir, "gse188646_pseudobulk_counts.mtx")
  rn_path <- file.path(out_dir, "gse188646_pseudobulk_counts_rownames.csv")
  cn_path <- file.path(out_dir, "gse188646_pseudobulk_counts_colnames.csv")
  meta_path <- file.path(out_dir, "gse188646_pseudobulk_metadata.csv")
  out_deg <- file.path(out_dir, "gse188646_young_vs_aged_deg.csv")
}

if (!file.exists(mtx_path)) stop("Missing mtx: ", mtx_path)
if (!file.exists(rn_path)) stop("Missing rownames csv: ", rn_path)
if (!file.exists(cn_path)) stop("Missing colnames csv: ", cn_path)
if (!file.exists(meta_path)) stop("Missing metadata csv: ", meta_path)

counts <- readMM(mtx_path)
rn <- read.csv(rn_path, stringsAsFactors = FALSE)
cn <- read.csv(cn_path, stringsAsFactors = FALSE)
gene_col <- if ("gene" %in% colnames(rn)) "gene" else colnames(rn)[[1]]
sample_col <- if ("sample" %in% colnames(cn)) "sample" else colnames(cn)[[1]]
rownames(counts) <- rn[[gene_col]]
colnames(counts) <- cn[[sample_col]]

meta <- read.csv(meta_path, stringsAsFactors = FALSE)
meta$age_bin <- factor(meta$age_bin, levels = c("Young", "Aged"))
meta <- meta[match(colnames(counts), meta$orig.ident), , drop = FALSE]
if (any(is.na(meta$orig.ident))) {
  stop("metadata orig.ident does not match count matrix columns")
}

y <- DGEList(counts = counts, group = meta$age_bin)
keep <- filterByExpr(y, group = meta$age_bin)
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

utils::write.csv(res, out_deg, row.names = FALSE)
message("Wrote DE from mtx: ", normalizePath(out_deg, winslash = "/", mustWork = FALSE))
