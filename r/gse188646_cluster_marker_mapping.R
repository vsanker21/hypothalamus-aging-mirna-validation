# Map seurat_clusters to putative cell classes using curated marker modules (mouse).
# No label transfer / reference atlas — descriptive mean expression per cluster.
#
# Usage (from feasibility_study/):
#   Rscript r/gse188646_cluster_marker_mapping.R path/to/GSE188646_hypo.integrated.final.20210719.RDS
#
# Outputs (outputs/gse188646_cluster_annotation/):
#   cluster_marker_means_long.csv   — cluster, gene, mean_expr
#   cluster_module_scores.csv       — cluster, module, mean_of_included_markers
#   cluster_putative_labels.csv     — cluster, rank1_module, rank1_z, rank2_module, rank2_z, delta_z
#   CLUSTER_MARKER_MAPPING_README.txt
#
# Caveats: integrated snRNA can distort absolute expression; modules overlap; rare populations
# may score low. Use alongside GEO figure annotations and orthogonal references (e.g. HypoMap).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

argv <- commandArgs(trailingOnly = FALSE)
fa <- grep("^--file=", argv, value = TRUE)
script_dir <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript r/gse188646_cluster_marker_mapping.R <seurat.rds>")
}
rds <- args[[1]]
if (!file.exists(rds)) stop("RDS not found: ", rds)

out_dir <- file.path(dirname(script_dir), "outputs", "gse188646_cluster_annotation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

obj <- readRDS(rds)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("Column seurat_clusters missing from metadata.")
}
if (!"RNA" %in% names(obj@assays)) {
  stop("RNA assay missing.")
}
DefaultAssay(obj) <- "RNA"

# Curated marker sets (MGI-style symbols; intersected with object at runtime).
# Sources: common hypothalamus / CNS snRNA atlases; illustrative not exhaustive.
modules <- list(
  Astrocyte = c("Gfap", "Slc1a2", "Aqp4", "Slc1a3", "Atp1b2"),
  Microglia = c("Cx3cr1", "P2ry12", "C1qa", "C1qb", "Aif1", "Hexb"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp"),
  OPC = c("Pdgfra", "Cspg4", "Sox10", "Vcan"),
  Endothelial = c("Pecam1", "Cldn5", "Vwf", "Ly6c1", "Kdr"),
  Pericyte = c("Pdgfrb", "Rgs5", "Acta2", "Cspg4"),
  Excitatory_neuron = c("Slc17a6", "Slc17a7", "Neurod1", "Satb2", "Tbr1"),
  Inhibitory_neuron = c("Gad1", "Gad2", "Slc32a1", "Slc6a1", "Pvalb", "Sst", "Vip"),
  Dopaminergic = c("Th", "Slc6a3", "Ddc", "Slc18a2"),
  Histaminergic = c("Hdc", "Slc22a3"),
  Serotonergic = c("Tph2", "Slc6a4", "Fev"),
  POMC = c("Pomc", "Pcsk1"),
  AGRP_NPY = c("Agrp", "Npy"),
  Kisspeptin = c("Kiss1", "Tac3"),
  Oxytocin = c("Oxt"),
  Vasopressin = c("Avp"),
  CRH = c("Crh", "Crhbp"),
  Tanycyte_ependymal = c("Foxj1", "Dcdc2b", "Rarres2", "Vim", "Rax", "Col23a1", "Ccdc153"),
  NSC_like = c("Sox2", "Nes", "Hes1", "Hes5", "Ascl1", "Prom1", "Fabp7"),
  Radial_glia_like = c("Vim", "Pax6", "Hes1", "Hes5", "Fabp7", "Slc1a3", "Gfap")
)

all_mk <- unique(unlist(modules))
rn <- rownames(obj)
feats <- intersect(all_mk, rn)
if (length(feats) < 10) {
  stop("Too few marker genes found in object (", length(feats), "). Check gene naming.")
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

pb <- avg_wrap(obj, feats, "seurat_clusters")
mat <- as.matrix(pb[["RNA"]])
if (is.null(colnames(mat))) {
  stop("AverageExpression returned no column names.")
}

# Long format: per-gene means
long <- data.frame(
  cluster = rep(colnames(mat), each = nrow(mat)),
  seurat_cluster_id = rep(sub("^g", "", colnames(mat)), each = nrow(mat)),
  gene = rep(rownames(mat), times = ncol(mat)),
  mean_expr = as.vector(mat),
  stringsAsFactors = FALSE
)
utils::write.csv(long, file.path(out_dir, "cluster_marker_means_long.csv"), row.names = FALSE)

# Module scores = row mean of member genes present in mat
mod_rows <- list()
for (nm in names(modules)) {
  g <- intersect(modules[[nm]], rownames(mat))
  if (length(g) == 0) {
    mod_rows[[nm]] <- rep(NA_real_, ncol(mat))
  } else {
    mod_rows[[nm]] <- colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
  }
}
mod_df <- as.data.frame(mod_rows, check.names = FALSE)
mod_df$cluster <- colnames(mat)
mod_df <- mod_df[, c("cluster", setdiff(names(mod_df), "cluster"))]
mod_df$seurat_cluster_id <- sub("^g", "", mod_df$cluster)
mod_df <- mod_df[, c("cluster", "seurat_cluster_id", setdiff(names(mod_df), c("cluster", "seurat_cluster_id")))]

# Z-score each module across clusters; rank1 / rank2 label
cn <- setdiff(colnames(mod_df), c("cluster", "seurat_cluster_id"))
zm <- as.matrix(mod_df[, cn, drop = FALSE])
for (j in seq_len(ncol(zm))) {
  v <- zm[, j]
  s <- stats::sd(v, na.rm = TRUE)
  zm[, j] <- if (!is.finite(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
}
lab <- data.frame(cluster = mod_df$cluster, stringsAsFactors = FALSE)
for (i in seq_len(nrow(zm))) {
  o <- order(zm[i, ], decreasing = TRUE, na.last = NA)
  lab$rank1_module[i] <- cn[o[[1]]]
  lab$rank1_z[i] <- zm[i, o[[1]]]
  if (length(o) >= 2) {
    lab$rank2_module[i] <- cn[o[[2]]]
    lab$rank2_z[i] <- zm[i, o[[2]]]
    lab$delta_z[i] <- zm[i, o[[1]]] - zm[i, o[[2]]]
  } else {
    lab$rank2_module[i] <- NA_character_
    lab$rank2_z[i] <- NA_real_
    lab$delta_z[i] <- NA_real_
  }
}
lab$seurat_cluster_id <- sub("^g", "", lab$cluster)
utils::write.csv(mod_df, file.path(out_dir, "cluster_module_scores.csv"), row.names = FALSE)

# Seurat may prefix cluster levels that start with digits (e.g. g0); align to numeric stratum IDs.
man_path <- file.path(dirname(script_dir), "outputs", "gse188646_strata", "manifest.csv")
if (file.exists(man_path)) {
  man <- utils::read.csv(man_path, stringsAsFactors = FALSE)
  if ("stratum" %in% colnames(man) && "stratum_col" %in% colnames(man)) {
    m2 <- man[man$stratum_col == "seurat_clusters", c("stratum", "n_cells", "n_young_rep", "n_aged_rep", "deg_csv")]
    colnames(m2)[1] <- "seurat_cluster_id"
    m2$seurat_cluster_id <- as.character(m2$seurat_cluster_id)
    lab$seurat_cluster_id <- as.character(lab$seurat_cluster_id)
    lab <- merge(lab, m2, by = "seurat_cluster_id", all.x = TRUE, sort = FALSE)
  }
}
tb <- table(obj$seurat_clusters)
tb <- tb[as.character(names(tb))]
if ("n_cells" %in% colnames(lab)) {
  miss <- is.na(lab$n_cells)
  if (any(miss)) {
    lab$n_cells[miss] <- as.numeric(tb[lab$seurat_cluster_id[miss]])
  }
} else {
  lab$n_cells <- as.numeric(tb[lab$seurat_cluster_id])
}
utils::write.csv(lab, file.path(out_dir, "cluster_putative_labels.csv"), row.names = FALSE)

readme <- c(
  "GSE188646 cluster → putative cell class (marker modules)",
  "=======================================================",
  "",
  "Method: Seurat::AverageExpression on normalized RNA (`data` layer), grouped by `seurat_clusters`.",
  "Each module score is the mean of detected marker genes in that cluster.",
  "Per-module z-scores are computed across clusters; rank1 = highest z (heuristic label).",
  "",
  "LIMITATIONS (read before interpreting DE strata):",
  "- Integrated object: batch/integration can shift relative expression between clusters.",
  "- Modules overlap (e.g. Pdgfra/Cspg4 appear in OPC and Pericyte lists); rank2 and delta_z indicate ambiguity.",
  "- Absence of a canonical `cell_type` column in GEO metadata — this is inference, not author labels.",
  "- Cross-dataset geometry: GSE208355 DE-axis correlations vs clusters are in outputs/gse188646_hypomap_mapping/ (separate script); not full atlas label transfer.",
  "",
  "Files:",
  "- cluster_marker_means_long.csv — raw mean normalized expression per gene × cluster",
  "- cluster_module_scores.csv — module means",
  "- cluster_putative_labels.csv — heuristic rank1/rank2 labels; seurat_cluster_id strips Seurat's leading-`g` prefix; n_cells from strata manifest when present else nucleus counts from full RDS table; deg_csv only for clusters that passed stratified DE filters",
  ""
)
writeLines(readme, file.path(out_dir, "CLUSTER_MARKER_MAPPING_README.txt"))

message("Wrote cluster annotation under: ", normalizePath(out_dir, winslash = "/"))
