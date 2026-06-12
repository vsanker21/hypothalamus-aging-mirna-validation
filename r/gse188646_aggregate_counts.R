# Shared: pseudobulk counts = sum RNA counts per orig.ident (chunked or one-shot).
# Chunked path partitions cells into blocks of size chunk_size to cap peak memory
# during Seurat::AggregateExpression on large snRNA-seq objects.
#
# Note: subassignment into acc must happen in this function's environment only (no nested
# closure mutating acc), otherwise R can leave acc all-zero and edgeR sees zero library sizes.

aggregate_pseudobulk_counts <- function(obj, chunk_size = 0L) {
  if (!"orig.ident" %in% colnames(obj@meta.data)) {
    stop("orig.ident required in metadata.")
  }
  genes <- rownames(obj)
  ocols <- sort(unique(as.character(obj$orig.ident)))
  cells <- colnames(obj)
  acc <- matrix(0, nrow = length(genes), ncol = length(ocols), dimnames = list(genes, ocols))

  run_agg <- function(o) {
    tryCatch(
      Seurat::AggregateExpression(
        o,
        assays = "RNA",
        slot = "counts",
        return.seurat = FALSE,
        group.by = "orig.ident",
        verbose = FALSE
      ),
      error = function(e) {
        Seurat::AggregateExpression(
          o,
          assays = "RNA",
          return.seurat = FALSE,
          group.by = "orig.ident",
          verbose = FALSE
        )
      }
    )
  }

  fix_cn <- function(cn) {
    gsub("^(Young|Aged)-([0-9]+)$", "\\1_\\2", cn, perl = TRUE)
  }

  if (is.null(chunk_size) || is.na(chunk_size) || chunk_size <= 0L) {
    pb <- run_agg(obj)
    mat <- as.matrix(pb[["RNA"]])
    rn <- rownames(mat)
    mm <- mat[match(genes, rn), , drop = FALSE]
    mm[is.na(mm)] <- 0
    cn <- fix_cn(colnames(mm))
    colnames(mm) <- cn
    for (j in colnames(mm)) {
      if (j %in% ocols) {
        acc[, j] <- acc[, j] + as.numeric(mm[, j])
      }
    }
  } else {
    n <- length(cells)
    message("Chunked AggregateExpression: ", n, " cells, chunk_size=", chunk_size)
    brk <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
    for (k in seq_along(brk)) {
      idx <- brk[[k]]
      sub <- subset(obj, cells = cells[idx])
      if (ncol(sub) == 0) next
      message("Chunk ", k, "/", length(brk), " (", ncol(sub), " cells) ...")
      pb <- run_agg(sub)
      mat <- as.matrix(pb[["RNA"]])
      rn <- rownames(mat)
      mm <- mat[match(genes, rn), , drop = FALSE]
      mm[is.na(mm)] <- 0
      cn <- fix_cn(colnames(mm))
      colnames(mm) <- cn
      for (j in colnames(mm)) {
        if (j %in% ocols) {
          acc[, j] <- acc[, j] + as.numeric(mm[, j])
        }
      }
      rm(sub, pb, mat, mm)
      gc(verbose = FALSE)
    }
  }

  cs <- colSums(acc)
  if (any(cs <= 0)) {
    bad <- names(cs)[cs <= 0]
    stop(
      "Zero total counts for sample(s): ",
      paste(bad, collapse = ", "),
      " — aggregation failed or metadata/colnames mismatch."
    )
  }

  Matrix::Matrix(acc, sparse = TRUE)
}
