# ==============================================================================
# integration.R - Multi-sample Integration Functions
# ==============================================================================
# Functions for integrating multiple scRNA-seq samples using:
# - Seurat reciprocal PCA (rPCA) for GBM cohort integration
# - BBKNN for batch correction (used in meta-atlas)
# ==============================================================================

#' Integrate samples using Seurat rPCA
#'
#' @param object_list List of Seurat objects
#' @param config Configuration list
#' @param split_by Metadata column for splitting (if single merged object)
#' @param verbose Print progress
#' @return Integrated Seurat object
#' @export
integrate_rpca <- function(object_list, config, split_by = NULL, 
                            verbose = TRUE) {
  
  if (verbose) log_message("Integrating samples using rPCA...")
  
  int_config <- config$integration
  
  # If single object, split by batch
  if (!is.list(object_list)) {
    if (is.null(split_by)) {
      stop("Must provide split_by column for single object")
    }
    object_list <- Seurat::SplitObject(object_list, split.by = split_by)
  }
  
  n_samples <- length(object_list)
  if (verbose) log_message(sprintf("  Processing %d samples", n_samples))
  
  # Normalize each sample if not already done
  object_list <- lapply(object_list, function(obj) {
    if (!"SCT" %in% names(obj@assays)) {
      obj <- Seurat::SCTransform(obj, verbose = FALSE)
    }
    return(obj)
  })
  
  # Select integration features
  features <- Seurat::SelectIntegrationFeatures(
    object_list,
    nfeatures = int_config$n_features
  )
  
  if (verbose) log_message(sprintf("  Using %d integration features", 
                                    length(features)))
  
  # Prep SCT integration
  object_list <- Seurat::PrepSCTIntegration(
    object_list,
    anchor.features = features
  )
  
  # Run PCA on each object
  object_list <- lapply(object_list, function(obj) {
    Seurat::RunPCA(obj, features = features, verbose = FALSE)
  })
  
  # Find integration anchors using rPCA
  anchors <- Seurat::FindIntegrationAnchors(
    object_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca",
    k.filter = int_config$k_filter,
    k.anchor = int_config$k_anchor,
    verbose = verbose
  )
  
  # Integrate data
  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    k.weight = int_config$k_weight,
    verbose = verbose
  )
  
  # Set default assay and run dimensionality reduction
  Seurat::DefaultAssay(integrated) <- "integrated"
  
  integrated <- Seurat::RunPCA(integrated, verbose = FALSE)
  n_pcs <- get_npcs_for_variance(integrated, 
                                  config$dim_reduction$pca$variance_threshold)
  
  integrated <- Seurat::RunUMAP(
    integrated,
    dims = 1:n_pcs,
    n.neighbors = config$dim_reduction$umap$n_neighbors,
    min.dist = config$dim_reduction$umap$min_dist,
    metric = config$dim_reduction$umap$metric,
    verbose = FALSE
  )
  
  integrated <- Seurat::FindNeighbors(integrated, dims = 1:n_pcs, verbose = FALSE)
  integrated <- Seurat::FindClusters(integrated, 
                                      resolution = config$clustering$resolution$coarse,
                                      verbose = FALSE)
  
  # Store n_pcs
  integrated@misc[["n_pcs"]] <- n_pcs
  
  if (verbose) {
    log_message(sprintf("  Integration complete: %d cells, %d clusters",
                        ncol(integrated),
                        length(unique(Seurat::Idents(integrated)))))
  }
  
  return(integrated)
}

#' Integrate using BBKNN (requires reticulate/Python)
#'
#' @param object Seurat object with multiple samples
#' @param batch_key Metadata column with batch info
#' @param config Configuration list
#' @param use_harmony Fallback to Harmony if BBKNN fails
#' @return Seurat object with corrected UMAP
#' @export
integrate_bbknn <- function(object, batch_key = "sample", config = NULL,
                             use_harmony = TRUE) {
  
  log_message("Integrating using BBKNN...")
  
  # Check if Python/BBKNN available
  bbknn_available <- tryCatch({
    reticulate::import("bbknn", convert = FALSE)
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  if (!bbknn_available) {
    if (use_harmony) {
      log_message("  BBKNN not available, falling back to Harmony...", "WARNING")
      return(integrate_harmony(object, batch_key, config))
    } else {
      stop("BBKNN not available and Harmony fallback disabled")
    }
  }
  
  # Import Python modules
  anndata <- reticulate::import("anndata", convert = FALSE)
  bbknn <- reticulate::import("bbknn", convert = FALSE)
  sc <- reticulate::import("scanpy", convert = FALSE)
  
  # Get PCA embeddings
  if (!"pca" %in% names(object@reductions)) {
    object <- Seurat::RunPCA(object, verbose = FALSE)
  }
  
  pca_embeddings <- Seurat::Embeddings(object, "pca")
  batch_labels <- object@meta.data[[batch_key]]
  
  # Create AnnData object
  adata <- anndata$AnnData(
    X = pca_embeddings,
    obs = data.frame(batch = batch_labels)
  )
  sc$tl$pca(adata)
  adata$obsm[["X_pca"]] <- pca_embeddings
  
  # Run BBKNN
  bbknn$bbknn(adata, batch_key = "batch")
  
  # Run UMAP on corrected neighbors
  sc$tl$umap(adata)
  
  # Extract UMAP coordinates
  umap_coords <- reticulate::py_to_r(adata$obsm[["X_umap"]])
  colnames(umap_coords) <- c("BBKNN_1", "BBKNN_2")
  rownames(umap_coords) <- colnames(object)
  
  # Add as new reduction
  object[["bbknn"]] <- Seurat::CreateDimReducObject(
    embeddings = umap_coords,
    key = "BBKNN_",
    assay = Seurat::DefaultAssay(object)
  )
  
  log_message("  BBKNN integration complete")
  
  return(object)
}

#' Integrate using Harmony as BBKNN fallback
#'
#' @param object Seurat object
#' @param batch_key Metadata column with batch info
#' @param config Configuration list
#' @return Seurat object with corrected embeddings
#' @export
integrate_harmony <- function(object, batch_key = "sample", config = NULL) {
  
  log_message("Integrating using Harmony...")
  
  # Run PCA if needed
  if (!"pca" %in% names(object@reductions)) {
    object <- Seurat::RunPCA(object, verbose = FALSE)
  }
  
  # Run Harmony
  object <- harmony::RunHarmony(
    object,
    group.by.vars = batch_key,
    reduction = "pca",
    assay.use = Seurat::DefaultAssay(object),
    verbose = FALSE
  )
  
  # Run UMAP on Harmony embeddings
  n_dims <- 30
  if (!is.null(config)) {
    n_dims <- config$integration$bbknn$n_pcs
  }
  
  object <- Seurat::RunUMAP(
    object,
    reduction = "harmony",
    dims = 1:n_dims,
    reduction.name = "umap_harmony",
    verbose = FALSE
  )
  
  log_message("  Harmony integration complete")
  
  return(object)
}

#' Merge multiple Seurat objects
#'
#' @param object_list List of Seurat objects
#' @param add_cell_ids Prefixes to add to cell names
#' @param project Project name
#' @return Merged Seurat object
#' @export
merge_seurat_objects <- function(object_list, add_cell_ids = NULL,
                                  project = "merged") {
  
  log_message(sprintf("Merging %d Seurat objects...", length(object_list)))
  
  if (is.null(add_cell_ids)) {
    add_cell_ids <- paste0("S", seq_along(object_list), "_")
  }
  
  # Merge
  merged <- object_list[[1]]
  
  if (length(object_list) > 1) {
    merged <- merge(
      merged,
      y = object_list[-1],
      add.cell.ids = add_cell_ids,
      project = project
    )
  }
  
  log_message(sprintf("  Merged object: %d cells, %d genes",
                      ncol(merged), nrow(merged)))
  
  return(merged)
}

#' Fix SCT for marker finding after integration
#'
#' @param object Seurat object after integration
#' @param assay Assay name (default "SCT")
#' @return Seurat object with fixed SCT
#' @export
fix_sct_for_markers <- function(object, assay = "SCT") {
  
  log_message("Preparing SCT for marker finding...")
  
  # Check if SCT models exist
  if (!assay %in% names(object@assays)) {
    warning("SCT assay not found")
    return(object)
  }
  
  # Get SCT models
  sct_models <- Seurat::SCTResults(object, slot = "cell.attributes")
  
  # Check for valid models
  valid_models <- vapply(sct_models, function(x) {
    !is.null(x) && nrow(x) >= 2
  }, logical(1))
  
  # Remove invalid models
  if (any(!valid_models)) {
    log_message(sprintf("  Removing %d invalid SCT models", sum(!valid_models)))
    object@assays[[assay]]@SCTModel.list <- 
      object@assays[[assay]]@SCTModel.list[valid_models]
  }
  
  # Prep for marker finding
  object <- Seurat::PrepSCTFindMarkers(object, assay = assay, verbose = FALSE)
  
  return(object)
}

#' Split integrated object by cell type for sub-clustering
#'
#' @param object Integrated Seurat object
#' @param cell_type_col Metadata column with cell type annotations
#' @param cell_types Cell types to extract (NULL for all)
#' @return Named list of Seurat objects
#' @export
split_by_celltype <- function(object, cell_type_col = "cell_type",
                               cell_types = NULL) {
  
  if (!cell_type_col %in% colnames(object@meta.data)) {
    stop(sprintf("Column '%s' not found in metadata", cell_type_col))
  }
  
  if (is.null(cell_types)) {
    cell_types <- unique(object@meta.data[[cell_type_col]])
  }
  
  result <- list()
  
  for (ct in cell_types) {
    cells <- colnames(object)[object@meta.data[[cell_type_col]] == ct]
    
    if (length(cells) < 10) {
      log_message(sprintf("  Skipping %s: only %d cells", ct, length(cells)),
                  "WARNING")
      next
    }
    
    result[[ct]] <- subset(object, cells = cells)
    log_message(sprintf("  %s: %d cells", ct, length(cells)))
  }
  
  return(result)
}

#' Re-embed subset of cells
#'
#' @param object Seurat object (subset)
#' @param config Configuration list
#' @param re_cluster Whether to re-cluster
#' @return Re-embedded Seurat object
#' @export
re_embed_subset <- function(object, config, re_cluster = TRUE) {
  
  log_message(sprintf("Re-embedding %d cells...", ncol(object)))
  
  # Set default assay
  if ("integrated" %in% names(object@assays)) {
    Seurat::DefaultAssay(object) <- "integrated"
  } else if ("SCT" %in% names(object@assays)) {
    Seurat::DefaultAssay(object) <- "SCT"
  }
  
  # Run PCA
  object <- Seurat::RunPCA(object, verbose = FALSE)
  n_pcs <- get_npcs_for_variance(object, 
                                  config$dim_reduction$pca$variance_threshold)
  
  # Run UMAP
  object <- Seurat::RunUMAP(
    object,
    dims = 1:n_pcs,
    n.neighbors = config$dim_reduction$umap$n_neighbors,
    min.dist = config$dim_reduction$umap$min_dist,
    metric = config$dim_reduction$umap$metric,
    verbose = FALSE
  )
  
  # Clustering
  if (re_cluster) {
    object <- Seurat::FindNeighbors(object, dims = 1:n_pcs, verbose = FALSE)
    object <- Seurat::FindClusters(
      object,
      resolution = config$clustering$resolution$fine,
      verbose = FALSE
    )
  }
  
  object@misc[["n_pcs"]] <- n_pcs
  
  return(object)
}
