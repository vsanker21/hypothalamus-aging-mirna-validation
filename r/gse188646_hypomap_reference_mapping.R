# Align GSE188646 seurat_clusters to HypoMap-related GSE208355 DESeq2 axes (GEO supplementary).
# NOT a download of the full HypoMap scRNA atlas (multi-GB); uses IP/input DE statistics as gene-level
# reference weights, plus optional user matrix (genes × cell types) via env GSE188646_HYPOMAP_REF_EXPR_CSV.
#
# Usage (from feasibility_study/):
#   Rscript r/gse188646_hypomap_reference_mapping.R <seurat.rds>
#
# Inputs:
#   - outputs/hypomap_geo_deseq_filtered_concat.csv (from run_extended.py) if present; else downloads
#     the three GSE208355 supplementary tables and applies the same padj / |log2FC| filters.
#   - Optional: Sys.setenv(GSE188646_HYPOMAP_REF_EXPR_CSV = "path/to/genes_x_celltypes.csv")
#     with column 'gene' (symbol) and one column per reference cell type (numeric means).
#
# Outputs (outputs/gse188646_hypomap_mapping/):
#   HYPomap_REFERENCE_MAPPING_README.txt
#   hypomap_axis_spearman.csv          — cluster × GSE208355 contrast (Spearman rho, p, n genes)
#   hypomap_axis_multiregress.csv      — per-cluster lm(expr ~ LFC_pomc + LFC_agrp + LFC_glp1r) across genes
#   hypomap_custom_ref_spearman.csv    — only if custom CSV present
#   cluster_hypomap_axes_combined.csv  — merge with cluster_putative_labels.csv when present

suppressPackageStartupMessages({
  library(Seurat)
})

argv <- commandArgs(trailingOnly = FALSE)
fa <- grep("^--file=", argv, value = TRUE)
script_dir <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
root <- dirname(script_dir)
out_dir <- file.path(root, "outputs", "gse188646_hypomap_mapping")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript r/gse188646_hypomap_reference_mapping.R <seurat.rds>")
}
rds <- args[[1]]
if (!file.exists(rds)) stop("RDS not found: ", rds)

padj_max <- 0.05
abs_lfc_min <- 0.5

urls <- list(
  pomc_deseq = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208355/suppl/GSE208355_pomc_deseq.csv.gz",
  agrp_deseq = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208355/suppl/GSE208355_agrp_deseq.csv.gz",
  glp1r_deseq = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208355/suppl/GSE208355_glp1r_deseq.csv.gz"
)

hypo_path <- file.path(root, "outputs", "hypomap_geo_deseq_filtered_concat.csv")
if (file.exists(hypo_path)) {
  message("Using existing ", hypo_path)
  hypo <- utils::read.csv(hypo_path, stringsAsFactors = FALSE)
} else {
  message("Building HypoMap GEO table from supplementary URLs (no local concat found)...")
  rows <- list()
  for (nm in names(urls)) {
    u <- urls[[nm]]
    tf <- tempfile(fileext = ".csv.gz")
    utils::download.file(u, tf, mode = "wb", quiet = TRUE)
    sub <- utils::read.csv(gzfile(tf), stringsAsFactors = FALSE)
    unlink(tf)
    sub$source_table <- nm
    if ("padj" %in% colnames(sub)) sub <- sub[!is.na(sub$padj) & sub$padj <= padj_max, , drop = FALSE]
    if ("log2FoldChange" %in% colnames(sub)) {
      sub <- sub[!is.na(sub$log2FoldChange) & abs(sub$log2FoldChange) >= abs_lfc_min, , drop = FALSE]
    }
    rows[[nm]] <- sub
  }
  hypo <- do.call(rbind, rows)
}

if (!"external_gene_name" %in% colnames(hypo)) {
  stop("HypoMap table missing external_gene_name.")
}
if (!"log2FoldChange" %in% colnames(hypo)) {
  stop("HypoMap table missing log2FoldChange.")
}
if (!"source_table" %in% colnames(hypo)) {
  hypo$source_table <- "combined"
}

# One weight per gene per contrast (mean log2FC if duplicated symbols)
weights_by_contrast <- list()
for (ctr in unique(hypo$source_table)) {
  h <- hypo[hypo$source_table == ctr, , drop = FALSE]
  g <- trimws(as.character(h$external_gene_name))
  w <- as.numeric(h$log2FoldChange)
  ok <- !is.na(g) & g != "" & !is.na(w)
  g <- g[ok]
  w <- w[ok]
  agg <- tapply(w, g, mean)
  weights_by_contrast[[ctr]] <- agg
}

obj <- readRDS(rds)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("seurat_clusters missing from RDS metadata.")
}
DefaultAssay(obj) <- "RNA"

# Genes: union of HypoMap DE genes present in RNA assay (keeps analysis focused & reproducible)
all_wgenes <- unique(unlist(lapply(weights_by_contrast, names)))
feats <- intersect(all_wgenes, rownames(obj))
if (length(feats) < 200) {
  stop("Too few intersecting genes (", length(feats), ") for stable mapping.")
}

avg_wrap <- function(o, features, group_by) {
  tryCatch(
    Seurat::AverageExpression(
      o,
      assays = "RNA",
      features = features,
      return.seurat = FALSE,
      group.by = group_by,
      layer = "data",
      verbose = FALSE
    ),
    error = function(e) {
      Seurat::AverageExpression(
        o,
        assays = "RNA",
        features = features,
        return.seurat = FALSE,
        group.by = group_by,
        slot = "data",
        verbose = FALSE
      )
    }
  )
}

message("AverageExpression on ", length(feats), " HypoMap-filtered genes × seurat_clusters ...")
pb <- avg_wrap(obj, feats, "seurat_clusters")
mat <- as.matrix(pb[["RNA"]])
colnames(mat) <- colnames(mat) # g0, ...
clusters <- colnames(mat)

spearman_rows <- list()
lm_rows <- list()

for (cl in clusters) {
  y <- mat[, cl, drop = TRUE]
  names(y) <- rownames(mat)
  # Design matrix: gene -> LFC in each contrast (0 if gene absent from that contrast)
  ctr_names <- names(weights_by_contrast)
  W <- matrix(0, nrow = length(feats), ncol = length(ctr_names), dimnames = list(feats, ctr_names))
  for (j in seq_along(ctr_names)) {
    wv <- weights_by_contrast[[ctr_names[j]]]
    W[feats, j] <- as.numeric(wv[feats])
    W[is.na(W[, j]), j] <- 0
  }
  yy <- y[feats]
  ok <- is.finite(yy)
  yy <- yy[ok]
  WW <- W[ok, , drop = FALSE]
  # Spearman per contrast: genes that appear in that contrast's filtered DE list only
  for (j in seq_along(ctr_names)) {
    gctr <- intersect(names(weights_by_contrast[[ctr_names[j]]]), names(yy))
    if (length(gctr) < 30) next
    yj <- yy[gctr]
    xj <- as.numeric(weights_by_contrast[[ctr_names[j]]][gctr])
    use <- is.finite(yj) & is.finite(xj)
    if (sum(use) < 30) next
    ct <- suppressWarnings(stats::cor.test(yj[use], xj[use], method = "spearman", exact = FALSE))
    spearman_rows[[length(spearman_rows) + 1L]] <- data.frame(
      cluster = cl,
      seurat_cluster_id = sub("^g", "", cl),
      contrast = ctr_names[j],
      rho = unname(ct$estimate),
      p = ct$p.value,
      n_genes = sum(use),
      stringsAsFactors = FALSE
    )
  }
  # Multivariate regression across genes (expr ~ LFC axes); coefficients are not "cell abundance"
  df <- data.frame(y = as.numeric(yy), WW, check.names = FALSE)
  if (nrow(df) > ncol(WW) + 5) {
    fit <- try(stats::lm(y ~ ., data = df), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      s <- summary(fit)
      f <- s$fstatistic
      p_lm <- if (!is.null(f) && length(f) >= 3) {
        stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
      } else {
        NA_real_
      }
      coefs <- stats::coef(fit)
      rd <- data.frame(
        cluster = cl,
        seurat_cluster_id = sub("^g", "", cl),
        n_genes = nrow(df),
        r_squared = s$r.squared,
        adj_r_squared = s$adj.r.squared,
        f_pvalue = as.numeric(p_lm),
        intercept = unname(coefs["(Intercept)"]),
        stringsAsFactors = FALSE
      )
      for (nm in ctr_names) {
        rd[[paste0("coef_", nm)]] <- if (nm %in% names(coefs)) unname(coefs[[nm]]) else NA_real_
      }
      lm_rows[[length(lm_rows) + 1L]] <- rd
    }
  }
}

sp <- do.call(rbind, spearman_rows)
if (is.null(sp) || nrow(sp) == 0) {
  stop("No Spearman results produced.")
}
sp$p_adj <- stats::p.adjust(sp$p, method = "BH")
utils::write.csv(sp, file.path(out_dir, "hypomap_axis_spearman.csv"), row.names = FALSE)

lm_df <- if (length(lm_rows)) do.call(rbind, lm_rows) else NULL
if (!is.null(lm_df)) {
  utils::write.csv(lm_df, file.path(out_dir, "hypomap_axis_multiregress.csv"), row.names = FALSE)
}

# Optional: user-supplied reference matrix (gene × cell types)
cref <- Sys.getenv("GSE188646_HYPOMAP_REF_EXPR_CSV", unset = "")
if (nzchar(cref) && file.exists(cref)) {
  message("Custom reference matrix: ", cref)
  ref <- utils::read.csv(cref, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene" %in% colnames(ref)) {
    warning("Custom ref CSV must have 'gene' column; skipping custom mapping.")
  } else {
    ref$gene <- trimws(as.character(ref$gene))
    tcols <- setdiff(colnames(ref), "gene")
    cust <- list()
    for (tc in tcols) {
      wv <- ref[[tc]]
      names(wv) <- ref$gene
      wv <- as.numeric(wv)
      okg <- !is.na(ref$gene) & !is.na(wv)
      agg <- tapply(wv[okg], ref$gene[okg], mean)
      cust[[tc]] <- agg
    }
    genes_c <- Reduce(intersect, lapply(cust, names))
    genes_c <- intersect(genes_c, rownames(obj))
    if (length(genes_c) < 200) {
      warning("Too few genes for custom ref; skipping.")
    } else {
      pb2 <- avg_wrap(obj, genes_c, "seurat_clusters")
      m2 <- as.matrix(pb2[["RNA"]])
      crrows <- list()
      for (cl in colnames(m2)) {
        y <- m2[, cl]
        names(y) <- rownames(m2)
        for (tc in names(cust)) {
          w <- cust[[tc]][genes_c]
          yy <- y[genes_c]
          use <- is.finite(yy) & is.finite(w)
          if (sum(use) < 30) next
          ct <- suppressWarnings(stats::cor.test(yy[use], w[use], method = "spearman", exact = FALSE))
          crrows[[length(crrows) + 1L]] <- data.frame(
            cluster = cl,
            seurat_cluster_id = sub("^g", "", cl),
            ref_celltype = tc,
            rho = unname(ct$estimate),
            p = ct$p.value,
            n_genes = sum(use),
            stringsAsFactors = FALSE
          )
        }
      }
      csp <- do.call(rbind, crrows)
      if (!is.null(csp) && nrow(csp)) {
        csp$p_adj <- stats::p.adjust(csp$p, method = "BH")
        utils::write.csv(csp, file.path(out_dir, "hypomap_custom_ref_spearman.csv"), row.names = FALSE)
      }
    }
  }
}

# Merge best GEO-axis hit with marker-based labels
lab_path <- file.path(root, "outputs", "gse188646_cluster_annotation", "cluster_putative_labels.csv")
if (file.exists(lab_path)) {
  lab <- utils::read.csv(lab_path, stringsAsFactors = FALSE)
  best <- sp[order(sp$cluster, sp$p_adj, -abs(sp$rho)), , drop = FALSE]
  best <- best[!duplicated(best$cluster), ]
  colnames(best)[colnames(best) == "contrast"] <- "hypomap_best_contrast"
  colnames(best)[colnames(best) == "rho"] <- "hypomap_best_rho"
  colnames(best)[colnames(best) == "p_adj"] <- "hypomap_best_padj"
  colnames(best)[colnames(best) == "n_genes"] <- "n_genes_hypomap_axis"
  m <- merge(lab, best[, c("cluster", "hypomap_best_contrast", "hypomap_best_rho", "hypomap_best_padj", "n_genes_hypomap_axis")],
    by = "cluster", all.x = TRUE
  )
  utils::write.csv(m, file.path(out_dir, "cluster_hypomap_axes_combined.csv"), row.names = FALSE)
}

readme <- c(
  "HypoMap-aligned reference mapping (GSE188646 clusters)",
  "========================================================",
  "",
  "WHAT THIS IS",
  "- Gene-level weights are log2FoldChange from GEO GSE208355 supplementary DESeq2 tables",
  "  (POMC, AGRP, GLP1R neuronal IP vs input contrasts — NOT whole-tissue young vs aged).",
  "- For each seurat_cluster, we average normalized RNA (`data` layer) across cells, restricted to",
  "  genes present in those filtered HypoMap DE lists and in the GSE188646 RNA assay.",
  "- Spearman correlation tests whether cluster mean expression covaries with each contrast axis.",
  "- With thousands of genes per test, p-values can be extreme even when |rho| is small; prioritise rho and biology, not p alone.",
  "- Multivariate linear model: per cluster, across genes, y = mean expr ~ LFC_pomc + LFC_agrp + LFC_glp1r",
  "  (descriptive geometry in gene space; not a generative causal model).",
  "",
  "WHAT THIS IS NOT",
  "- Not label transfer from the full HypoMap scRNA atlas (multi-GB AnnData / cellxgene object).",
  "- Not cell-type identity: significant correlation with e.g. pomc_deseq does not imply 'POMC neurons'.",
  "- Cross-study differences (sex, age, sorting, batch, integration) apply — treat as hypothesis generation.",
  "",
  "OPTIONAL FULL ATLAS PROFILES",
  "- Set GSE188646_HYPOMAP_REF_EXPR_CSV to a CSV with column `gene` and one column per reference cell type",
  "  (mean log-normalized expression). If present, Spearman cluster-vs-type correlations are written to",
  "  hypomap_custom_ref_spearman.csv. You must build/supply that matrix from HypoMap or your own summary.",
  "",
  "FILES",
  "- hypomap_axis_spearman.csv",
  "- hypomap_axis_multiregress.csv (if lm succeeded per cluster)",
  "- hypomap_custom_ref_spearman.csv (optional)",
  "- cluster_hypomap_axes_combined.csv (merge with marker-module labels when available)",
  ""
)
writeLines(readme, file.path(out_dir, "HYPomap_REFERENCE_MAPPING_README.txt"))

message("Wrote HypoMap-axis mapping under: ", normalizePath(out_dir, winslash = "/"))
