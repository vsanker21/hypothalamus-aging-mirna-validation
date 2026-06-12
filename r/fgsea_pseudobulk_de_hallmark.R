# fgsea / fgseaMultilevel on pseudobulk edgeR DE table (gene-level ranks), same Hallmark GMT.
# Prefers fgseaMultilevel (adaptive multilevel sampling); falls back to fgsea on error.
# Exploratory: ranks mix signal and noise; not interchangeable with miRNA-target ORA.
#
# Usage (from feasibility_study/):
#   Rscript r/fgsea_pseudobulk_de_hallmark.R <deg_csv> <hallmark.gmt> <out_csv> [nperm]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript r/fgsea_pseudobulk_de_hallmark.R <deg_csv> <hallmark.gmt> <out_csv> [nperm]")
}
deg_csv <- args[[1]]
gmt_path <- args[[2]]
out_csv <- args[[3]]
nperm <- if (length(args) >= 4) as.integer(args[[4]]) else 2000L

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

df <- data.table::fread(deg_csv, data.table = FALSE)
if (!all(c("gene", "logFC") %in% colnames(df))) {
  stop("deg_csv must contain columns: gene, logFC")
}
df$gene <- toupper(as.character(df$gene))
pc <- if ("p_val" %in% colnames(df)) "p_val" else if ("PValue" %in% colnames(df)) "PValue" else NULL
if (is.null(pc)) stop("deg_csv must contain p_val or PValue")
pv <- pmax(as.numeric(df[[pc]]), 1e-300)
sc <- sign(as.numeric(df$logFC)) * (-log10(pv))
names(sc) <- as.character(df$gene)
sc <- sort(sc, decreasing = TRUE)
sc <- sc[!duplicated(names(sc))]

pathways <- fgsea::gmtPathways(gmt_path)
pathways <- lapply(pathways, function(g) toupper(as.character(g)))
# Restrict ranks to genes present in the GMT (reduces vector length; avoids fgsea native
# instability on very long rank vectors on some Windows R builds). Full-genome ranks:
#   Sys.setenv(GSE188646_FGSEA_FULL_RANK = "1") before Rscript.
use_full_rank <- Sys.getenv("GSE188646_FGSEA_FULL_RANK", unset = "") != ""
if (!use_full_rank) {
  ug <- unique(unlist(pathways, use.names = FALSE))
  ug <- intersect(ug, names(sc))
  if (length(ug) < 100L) {
    stop("Too few DE genes overlap Hallmark GMT symbols; check species / symbol matching.")
  }
  sc <- sc[names(sc) %in% ug]
  sc <- sort(sc, decreasing = TRUE)
}

set.seed(43)
min_sz <- 10L
max_sz <- 800L
out <- tryCatch(
  {
    r <- fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = sc,
      sampleSize = 101L,
      minSize = min_sz,
      maxSize = max_sz,
      scoreType = "std",
      nPermSimple = nperm,
      nproc = 1L
    )
    list(df = r, meth = "fgseaMultilevel")
  },
  error = function(e) {
    message("fgseaMultilevel failed; falling back to fgsea: ", conditionMessage(e))
    r2 <- fgsea::fgsea(
      pathways = pathways,
      stats = sc,
      minSize = min_sz,
      maxSize = max_sz,
      nPermSimple = nperm,
      scoreType = "std"
    )
    list(df = r2, meth = "fgsea")
  }
)
res <- as.data.frame(out$df)
if (nrow(res) == 0L) {
  message("fgseaMultilevel returned 0 rows; falling back to fgsea")
  res <- as.data.frame(
    fgsea::fgsea(
      pathways = pathways,
      stats = sc,
      minSize = min_sz,
      maxSize = max_sz,
      nPermSimple = nperm,
      scoreType = "std"
    )
  )
  meth_used <- "fgsea"
} else {
  meth_used <- out$meth
}
if (nrow(res) > 0L) {
  res$method <- meth_used
}
if ("leadingEdge" %in% colnames(res)) {
  res$leadingEdge <- vapply(
    res$leadingEdge,
    function(x) paste(as.character(unlist(x)), collapse = ";"),
    character(1)
  )
}
utils::write.csv(res, out_csv, row.names = FALSE)
meth <- if ("method" %in% colnames(res) && nrow(res) > 0L) unique(as.character(res$method)) else meth_used
message("Wrote ", meth[1], " (pseudobulk DE rank): ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
