# Cohort 2 (contrast-matched public anchor): GSE87102 C57BL/6 whole hypothalamus microarray,
# young (2–3 mo) vs aged (22–24 mo), male animals (GEO titles). Harmonised to gene symbols.
# Complements GSE188646 female snRNA pseudobulk (different sex/assay — disclosed in meta README).
#
# Requires: BiocManager packages GEOquery, limma
#
# Usage (from feasibility_study/):
#   Rscript r/cohort2_gse87102_c57_hypothalamus_limma.R
#
# Output:
#   outputs/cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv
#     columns: gene, logFC, se_logFC, p_val, p_val_adj (Aged vs Young, log2 microarray)

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery", ask = FALSE, update = FALSE, quiet = TRUE)
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma", ask = FALSE, update = FALSE, quiet = TRUE)
  }
  library(GEOquery)
  library(limma)
})

out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_csv <- file.path(out_dir, "cohort2_GSE87102_C57_hypothalamus_aged_vs_young_limma.csv")

es <- getGEO("GSE87102", GSEMatrix = TRUE)[[1]]
titles <- as.character(pData(es)$title)
keep <- grepl("^C57BL6-hypothalamus-(young|old)-", titles, perl = TRUE) &
  !grepl("middle", titles, ignore.case = TRUE)
es <- es[, keep]
titles <- as.character(pData(es)$title)
age <- ifelse(grepl("-young-", titles, fixed = TRUE), "Young", "Aged")
if (!all(age %in% c("Young", "Aged"))) {
  stop("Unexpected sample titles after filter.")
}

sym <- as.character(fData(es)$GENE_SYMBOL)
ok <- !is.na(sym) & sym != "" & sym != "---"
es <- es[ok, ]
sym <- as.character(fData(es)$GENE_SYMBOL)
ex <- limma::avereps(exprs(es), ID = sym)

group <- factor(age, levels = c("Young", "Aged"))
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Aged")

fit <- lmFit(ex, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "Aged", number = Inf, sort.by = "none")
se_lfc <- abs(tt$logFC) / pmax(abs(tt$t), .Machine$double.eps)
res <- data.frame(
  gene = rownames(tt),
  logFC = tt$logFC,
  se_logFC = as.numeric(se_lfc),
  p_val = tt$P.Value,
  p_val_adj = tt$adj.P.Val,
  stringsAsFactors = FALSE
)

utils::write.csv(res, out_csv, row.names = FALSE)
message("Wrote cohort2 DE: ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
