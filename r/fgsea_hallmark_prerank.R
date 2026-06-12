# Preranked fgsea vs MSigDB Hallmark GMT (same file as Python gseapy path).
# Usage (from feasibility_study/):
#   Rscript r/fgsea_hallmark_prerank.R <ranked_csv> <hallmark.gmt> <out_csv> [nperm]
#
# ranked_csv must contain columns: gene, rank_metric (higher = more ranked toward top).
#
# Requires Bioconductor package fgsea (installed on first run via BiocManager).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript r/fgsea_hallmark_prerank.R <ranked_csv> <hallmark.gmt> <out_csv> [nperm]")
}
rank_csv <- args[[1]]
gmt_path <- args[[2]]
out_csv <- args[[3]]
nperm <- if (length(args) >= 4) as.integer(args[[4]]) else 2000L
if (is.na(nperm) || nperm < 100L) nperm <- 1000L

if (!file.exists(rank_csv)) stop("ranked_csv not found: ", rank_csv)
if (!file.exists(gmt_path)) stop("gmt not found: ", gmt_path)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE, update = FALSE, quiet = TRUE)
}
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

rnk <- data.table::fread(rank_csv, data.table = FALSE)
if (!all(c("gene", "rank_metric") %in% colnames(rnk))) {
  stop("ranked_csv must contain columns: gene, rank_metric")
}
rnk <- rnk[!is.na(rnk$gene) & !is.na(rnk$rank_metric), , drop = FALSE]
rnk$gene <- as.character(rnk$gene)
# fgsea expects decreasing stats; keep rank_metric as provided (Python uses higher = stronger)
ord <- order(rnk$rank_metric, decreasing = TRUE, method = "radix")
rnk <- rnk[ord, , drop = FALSE]
dup <- duplicated(rnk$gene)
if (any(dup)) {
  rnk <- rnk[!dup, , drop = FALSE]
}
stats <- rnk$rank_metric
names(stats) <- rnk$gene

pathways <- fgsea::gmtPathways(gmt_path)
set.seed(42)
res <- fgsea::fgsea(
  pathways = pathways,
  stats = stats,
  minSize = 10L,
  maxSize = 800L,
  nPermSimple = nperm,
  scoreType = "pos"
)
res <- as.data.frame(res)
if ("leadingEdge" %in% colnames(res)) {
  res$leadingEdge <- vapply(
    res$leadingEdge,
    function(x) paste(as.character(unlist(x)), collapse = ";"),
    character(1)
  )
}
utils::write.csv(res, out_csv, row.names = FALSE)
message("Wrote fgsea results: ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
