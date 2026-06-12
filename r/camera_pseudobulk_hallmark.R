# limma::camera on pseudobulk log-CPM (same samples as edgeR pseudobulk), MSigDB Hallmark GMT.
# Exploratory gene-set test conditioned on the experimental design; complements fgsea on ranks.
#
# Usage (from feasibility_study/):
#   Rscript r/camera_pseudobulk_hallmark.R <counts.mtx> <rownames.csv> <colnames.csv> <metadata.csv> <hallmark.gmt> <out_csv>
#
# Requires: Matrix, edgeR, limma, fgsea (for gmtPathways)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(
    "Usage: Rscript r/camera_pseudobulk_hallmark.R ",
    "<counts.mtx> <rownames.csv> <colnames.csv> <metadata.csv> <hallmark.gmt> <out_csv>"
  )
}
mtx_path <- args[[1]]
rn_path <- args[[2]]
cn_path <- args[[3]]
meta_path <- args[[4]]
gmt_path <- args[[5]]
out_csv <- args[[6]]

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(fgsea)
})

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

y <- DGEList(counts = counts)
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 1)
design <- model.matrix(~ age_bin, data = meta)
colnames(design) <- sub("^age_bin", "", colnames(design))

pathways <- gmtPathways(gmt_path)
idx <- ids2indices(pathways, rownames(logcpm), remove.empty = TRUE)
cam <- camera(logcpm, idx, design, contrast = ncol(design))
out <- as.data.frame(cam)
out$pathway <- rownames(out)
utils::write.csv(out, out_csv, row.names = FALSE)
message("Wrote camera: ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
