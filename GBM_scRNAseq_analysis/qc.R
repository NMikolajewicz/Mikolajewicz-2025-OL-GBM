# ==============================================================================
# qc.R - Quality Control and Preprocessing Functions
# ==============================================================================
# Functions for scRNA-seq quality control, filtering, and normalization.
# 
# QC thresholds:
# - 200-9000 genes/cell
# - <60% mitochondrial content
# - 3-MAD unmatched-rate filter
# ==============================================================================

#' Calculate QC metrics for Seurat object
#'
#' @param object Seurat object
#' @param mito_pattern Pattern for mitochondrial genes
#' @param ribo_pattern Pattern for ribosomal genes
#' @return Seurat object with QC metrics in metadata
#' @export
calculate_qc_metrics <- function(object, 
                                  mito_pattern = "^MT-",
                                  ribo_pattern = "^RP[SL]") {
  
  log_message("Calculating QC metrics...")
  
  # Mitochondrial content
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object, pattern = mito_pattern
  )

  # Ribosomal content
  object[["percent.ribo"]] <- Seurat::PercentageFeatureSet(
    object, pattern = ribo_pattern
  )
  
  # Log10 UMI counts
  object[["log10_nCount"]] <- log10(object$nCount_RNA + 1)
  
  # Log10 gene counts
  object[["log10_nFeature"]] <- log10(object$nFeature_RNA + 1)
  
  # Genes per UMI (complexity)
  object[["complexity"]] <- object$nFeature_RNA / object$nCount_RNA
  
  log_message(sprintf("  Median genes/cell: %d", 
                      median(object$nFeature_RNA)))
  log_message(sprintf("  Median UMIs/cell: %d", 
                      median(object$nCount_RNA)))
  log_message(sprintf("  Median mito%%: %.1f", 
                      median(object$percent.mt)))
  
  return(object)
}

#' Apply QC filters based on configuration
#'
#' @param object Seurat object with QC metrics
#' @param config Configuration list with QC parameters
#' @param verbose Print filtering statistics
#' @return Filtered Seurat object
#' @export
apply_qc_filters <- function(object, config, verbose = TRUE) {
  
  n_cells_before <- ncol(object)
  
  # Get thresholds from config
  min_genes <- config$qc$min_genes
  max_genes <- config$qc$max_genes
  max_mito <- config$qc$max_mito_percent
  
  if (verbose) {
    log_message("Applying QC filters:")
    log_message(sprintf("  min_genes: %d", min_genes))
    log_message(sprintf("  max_genes: %d", max_genes))
    log_message(sprintf("  max_mito%%: %d", max_mito))
  }
  
  # Apply filters
  cells_pass <- object$nFeature_RNA >= min_genes &
                object$nFeature_RNA <= max_genes &
                object$percent.mt < max_mito
  
  object <- subset(object, cells = colnames(object)[cells_pass])
  
  n_cells_after <- ncol(object)
  
  if (verbose) {
    log_message(sprintf("  Cells before: %d", n_cells_before))
    log_message(sprintf("  Cells after: %d", n_cells_after))
    log_message(sprintf("  Cells removed: %d (%.1f%%)", 
                        n_cells_before - n_cells_after,
                        (n_cells_before - n_cells_after) / n_cells_before * 100))
  }
  
  return(object)
}

#' Apply MAD-based outlier detection (for doublet filtering)
#'
#' @param object Seurat object
#' @param metric Metric to use for MAD filtering
#' @param n_mad Number of MADs for threshold
#' @param direction "both", "upper", or "lower"
#' @return Logical vector of cells passing filter
#' @export
mad_filter <- function(object, metric, n_mad = 3, direction = "both") {
  
  values <- object@meta.data[[metric]]
  median_val <- median(values, na.rm = TRUE)
  mad_val <- mad(values, na.rm = TRUE)
  
  if (direction == "both") {
    pass <- values >= (median_val - n_mad * mad_val) &
            values <= (median_val + n_mad * mad_val)
  } else if (direction == "upper") {
    pass <- values <= (median_val + n_mad * mad_val)
  } else if (direction == "lower") {
    pass <- values >= (median_val - n_mad * mad_val)
  }
  
  return(pass)
}

#' Filter genes by expression
#'
#' @param object Seurat object
#' @param min_cells Minimum cells expressing gene
#' @param min_pct Minimum percentage of cells expressing gene
#' @return Vector of genes passing filter
#' @export
filter_genes <- function(object, min_cells = 3, min_pct = NULL) {
  
  n_cells <- ncol(object)
  
  # Calculate cells per gene
  counts <- Seurat::GetAssayData(object, slot = "counts")
  cells_per_gene <- Matrix::rowSums(counts > 0)
  
  if (!is.null(min_pct)) {
    min_cells <- ceiling(n_cells * min_pct)
  }
  
  genes_pass <- names(cells_per_gene)[cells_per_gene >= min_cells]
  
  log_message(sprintf("Genes passing filter (>=%d cells): %d/%d",
                      min_cells, length(genes_pass), nrow(object)))
  
  return(genes_pass)
}

#' SCTransform normalization with manuscript parameters
#'
#' @param object Seurat object
#' @param config Configuration list
#' @param vars_to_regress Variables to regress out (NULL to use config)
#' @param verbose Print progress
#' @return Normalized Seurat object
#' @export
normalize_sctransform <- function(object, config, 
                                   vars_to_regress = NULL,
                                   verbose = TRUE) {
  
  if (verbose) log_message("Running SCTransform normalization...")
  
  # Get parameters from config
  sct_config <- config$normalization$sctransform
  
  if (is.null(vars_to_regress)) {
    vars_to_regress <- sct_config$vars_to_regress
  }
  
  # Check if percent.mt exists if we want to regress it
  if ("percent.mt" %in% vars_to_regress && 
      !"percent.mt" %in% colnames(object@meta.data)) {
    log_message("  Adding percent.mt for regression...")
    object <- calculate_qc_metrics(object)
  }
  
  # Run SCTransform
  object <- tryCatch({
    Seurat::SCTransform(
      object,
      method = sct_config$method,
      vst.flavor = sct_config$vst_flavor,
      vars.to.regress = vars_to_regress,
      variable.features.n = sct_config$variable_features_n,
      verbose = verbose,
      conserve.memory = sct_config$conserve_memory
    )
  }, error = function(e) {
    log_message("  vst.flavor v2 failed, falling back to default...", "WARNING")
    Seurat::SCTransform(
      object,
      method = sct_config$method,
      vars.to.regress = vars_to_regress,
      variable.features.n = sct_config$variable_features_n,
      verbose = verbose
    )
  })
  
  if (verbose) {
    log_message(sprintf("  Variable features: %d", 
                        length(Seurat::VariableFeatures(object))))
  }
  
  return(object)
}

#' Run PCA with variance-based dimensionality selection
#'
#' @param object Seurat object (normalized)
#' @param config Configuration list
#' @param verbose Print progress
#' @return Seurat object with PCA
#' @export
run_pca <- function(object, config, verbose = TRUE) {
  
  if (verbose) log_message("Running PCA...")
  
  npcs <- config$dim_reduction$pca$npcs
  var_threshold <- config$dim_reduction$pca$variance_threshold
  
  # Run PCA
  object <- Seurat::RunPCA(object, npcs = npcs, verbose = FALSE)
  
  # Determine optimal PCs
  n_pcs <- get_npcs_for_variance(object, var_threshold)
  
  # Store in metadata
  object@misc[["n_pcs"]] <- n_pcs
  
  return(object)
}

#' Run UMAP with manuscript parameters
#'
#' @param object Seurat object with PCA
#' @param config Configuration list
#' @param dims PCA dimensions to use (NULL to use stored n_pcs)
#' @param verbose Print progress
#' @return Seurat object with UMAP
#' @export
run_umap <- function(object, config, dims = NULL, verbose = TRUE) {
  
  if (verbose) log_message("Running UMAP...")
  
  umap_config <- config$dim_reduction$umap
  
  # Get dimensions
  if (is.null(dims)) {
    if (!is.null(object@misc[["n_pcs"]])) {
      n_pcs <- object@misc[["n_pcs"]]
    } else {
      n_pcs <- get_npcs_for_variance(object, 
                                      config$dim_reduction$pca$variance_threshold)
    }
    dims <- 1:n_pcs
  }
  
  # Run UMAP with manuscript parameters
  object <- Seurat::RunUMAP(
    object,
    dims = dims,
    n.neighbors = umap_config$n_neighbors,
    min.dist = umap_config$min_dist,
    metric = umap_config$metric,
    verbose = FALSE
  )
  
  if (verbose) {
    log_message(sprintf("  Using %d PCA dimensions", length(dims)))
    log_message(sprintf("  n_neighbors: %d, min_dist: %.2f, metric: %s",
                        umap_config$n_neighbors,
                        umap_config$min_dist,
                        umap_config$metric))
  }
  
  return(object)
}

#' Run clustering with specified resolution
#'
#' @param object Seurat object with PCA
#' @param config Configuration list
#' @param resolution Clustering resolution (NULL to use config$coarse)
#' @param dims PCA dimensions to use
#' @param verbose Print progress
#' @return Seurat object with clusters
#' @export
run_clustering <- function(object, config, resolution = NULL,
                            dims = NULL, verbose = TRUE) {
  
  if (verbose) log_message("Running clustering...")
  
  # Get dimensions
  if (is.null(dims)) {
    if (!is.null(object@misc[["n_pcs"]])) {
      n_pcs <- object@misc[["n_pcs"]]
    } else {
      n_pcs <- 30
    }
    dims <- 1:n_pcs
  }
  
  # Get resolution
  if (is.null(resolution)) {
    resolution <- config$clustering$resolution$coarse
  }
  
  # Find neighbors
  object <- Seurat::FindNeighbors(
    object, 
    dims = dims,
    k.param = config$clustering$k_neighbors,
    verbose = FALSE
  )
  
  # Find clusters
  object <- Seurat::FindClusters(
    object,
    resolution = resolution,
    verbose = FALSE
  )
  
  n_clusters <- length(unique(Seurat::Idents(object)))
  
  if (verbose) {
    log_message(sprintf("  Resolution: %.2f", resolution))
    log_message(sprintf("  Number of clusters: %d", n_clusters))
  }
  
  return(object)
}

#' Complete preprocessing pipeline
#'
#' @param object Seurat object (raw counts)
#' @param config Configuration list
#' @param sample_name Sample identifier for logging
#' @return Preprocessed Seurat object
#' @export
preprocess_sample <- function(object, config, sample_name = "sample") {
  
  log_message(sprintf("\n=== Preprocessing: %s ===", sample_name))
  log_message(sprintf("Input cells: %d, genes: %d", ncol(object), nrow(object)))
  
  # QC metrics
  object <- calculate_qc_metrics(object)
  
  # QC filtering
  object <- apply_qc_filters(object, config)
  
  # Gene filtering
  genes_pass <- filter_genes(object, 
                              min_cells = config$qc$min_cells_per_gene)
  object <- subset(object, features = genes_pass)
  
  # Normalization
  object <- normalize_sctransform(object, config)
  
  # PCA
  object <- run_pca(object, config)
  
  # UMAP
  object <- run_umap(object, config)
  
  # Clustering
  object <- run_clustering(object, config)
  
  log_message(sprintf("Output cells: %d, genes: %d\n", ncol(object), nrow(object)))
  
  return(object)
}

#' Get expressed genes above threshold
#'
#' @param object Seurat object
#' @param min_pct Minimum percentage of cells expressing gene
#' @return Vector of gene names
#' @export
get_expressed_genes <- function(object, min_pct = 0.05) {
  counts <- Seurat::GetAssayData(object, slot = "counts")
  pct_expressing <- Matrix::rowMeans(counts > 0)
  genes <- names(pct_expressing)[pct_expressing >= min_pct]
  return(genes)
}

#' Plot QC metrics
#'
#' @param object Seurat object with QC metrics
#' @param config Configuration list
#' @return ggplot object
#' @export
plot_qc_metrics <- function(object, config = NULL) {
  
  df <- object@meta.data
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = nFeature_RNA)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::labs(x = "Genes per cell", y = "Count", title = "Gene Distribution") +
    ggplot2::theme_minimal()
  
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = nCount_RNA)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "UMIs per cell (log10)", y = "Count", title = "UMI Distribution") +
    ggplot2::theme_minimal()
  
  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = percent.mt)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::labs(x = "Mitochondrial %", y = "Count", title = "Mito Distribution") +
    ggplot2::theme_minimal()
  
  p4 <- ggplot2::ggplot(df, ggplot2::aes(x = nFeature_RNA, y = percent.mt)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.5) +
    ggplot2::labs(x = "Genes per cell", y = "Mitochondrial %", 
                  title = "Genes vs Mito") +
    ggplot2::theme_minimal()
  
  # Add threshold lines if config provided
  if (!is.null(config)) {
    p1 <- p1 + 
      ggplot2::geom_vline(xintercept = c(config$qc$min_genes, config$qc$max_genes),
                          linetype = "dashed", color = "red")
    p3 <- p3 +
      ggplot2::geom_vline(xintercept = config$qc$max_mito_percent,
                          linetype = "dashed", color = "red")
  }
  
  cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)
}
