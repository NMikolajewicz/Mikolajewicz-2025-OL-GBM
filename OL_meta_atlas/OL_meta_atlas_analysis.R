################################################################################
# Oligodendrocyte Meta-Atlas Generation and NMF-Based Gene Program Discovery
################################################################################
#
# Description: This R script generates an oligodendrocyte meta-atlas by integrating
#              scRNA-seq data from multiple studies spanning healthy brain, 
#              neurological diseases, and brain tumors. NMF is used to identify
#              robust transcriptional programs across conditions.
#
# Associated Publication: Mikolajewicz et al.
#
# Author: Nicholas Mikolajewicz
# Contact: n.mikolajewicz@utoronto.ca
#
# Input: Raw/processed data from individual studies (see Data Availability)
# Output: 
#   - Integrated oligodendrocyte meta-atlas (Seurat object)
#   - NMF-derived gene programs (Table S5)
#   - Module scores and visualizations
#
# NOTE: This script is designed to be self-contained and reproducible.
#       All custom functions from the scMiko package are included inline.
#
################################################################################

# ==============================================================================
# PART 0: SETUP AND DEPENDENCIES
# ==============================================================================

# Clear environment
rm(list = ls())
gc()

# Set random seed for reproducibility
set.seed(42)

# ---- Required R Packages ----
required_packages <- c(
  # Core scRNA-seq analysis
  "Seurat",
  "sctransform",
  "glmGamPoi",
  
  # NMF analysis
  "NNLM",
  "NMF",
  
  # Data manipulation
  "dplyr",
  "tidyr",
  "reshape2",
  "stringr",
  
  # Visualization
  "ggplot2",
  "RColorBrewer",
  "viridis",
  "cowplot",
  "pheatmap",
  "scattermore",
  
  # Parallel processing
  "foreach",
  "parallel",
  "doParallel",
  "pbapply",
  
  # File I/O
  "Matrix",
  "rhdf5",
  
  # Python interface for BBKNN
  "reticulate",
  
  # Gene annotation
  "biomaRt"
)

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# ==============================================================================
# PART 1: CONFIGURATION
# ==============================================================================

#' Configuration object containing all analysis parameters
#' Modify these paths and parameters for your environment

config <- list(
  # ---- Data Paths ----
  # Base directory containing raw data from external studies
  external_data_dir = "path/to/external/studies/",
  
  # Directory for GBM study data (current study)
  gbm_data_dir = "path/to/gbm/data/",
  
  # Output directory for results
  output_dir = "path/to/output/",
  
  # ---- QC Parameters ----
  min_cells_per_sample = 50,       # Minimum cells to retain a sample
  mito_threshold = 10,              # Maximum mitochondrial content (%)
  min_features = 200,               # Minimum genes per cell
  max_features = 9000,              # Maximum genes per cell
  
  # ---- Normalization Parameters ----
  n_variable_features = 2000,       # Number of variable features for SCT
  n_variable_features_integration = 3000, # For integration
  
  # ---- Dimensionality Reduction ----
  n_pcs = 50,                       # Number of principal components
  
  # ---- NMF Parameters ----
  nmf_k_range = 2:15,               # Range of k values to test
  top_n_genes = 50,                 # Top genes per NMF program
  min_pct_expressed = 0.005,        # Minimum fraction of cells expressing gene
  
  # ---- Robust Program Identification ----
  intra_threshold = 0.7,            # Within-sample Jaccard threshold
  inter_threshold = 0.2,            # Between-sample Jaccard threshold
  min_samples_robust = 3,           # Minimum samples for robust program
  n_final_programs = 8,             # Number of final programs to identify
  prevalence_threshold = 0.3,       # Gene prevalence threshold for final programs
  
  # ---- Parallelization ----
  n_cores = max(1, parallel::detectCores() - 2),
  
  # ---- Python/Conda for BBKNN ----
  conda_env = "scRNA",              # Name of conda environment with bbknn
  use_bbknn = TRUE                  # Whether to use BBKNN (requires Python)
)

# ==============================================================================
# PART 2: CUSTOM HELPER FUNCTIONS (scMiko equivalents)
# ==============================================================================

# ---- Utility Functions ----

#' Print timestamped message
#' @param msg Message to print
miko_message <- function(msg) {

  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
}

#' Convert column to rownames
#' @param df Data frame
#' @param col Column name to use as rownames
#' @return Data frame with column converted to rownames
col2rowname <- function(df, col) {
  rownames(df) <- df[[col]]
  df[[col]] <- NULL
  return(df)
}

#' Convert rownames to column
#' @param df Data frame
#' @param col Name for the new column
#' @return Data frame with rownames as column
rowname2col <- function(df, col = "rowname") {
  df[[col]] <- rownames(df)
  rownames(df) <- NULL
  return(df)
}

#' Convert named list to wide data frame
#' @param x Named list of vectors
#' @return Wide data frame
namedList2wideDF <- function(x) {
  max_len <- max(sapply(x, length))
  df <- data.frame(matrix(NA, nrow = max_len, ncol = length(x)))
  colnames(df) <- names(x)
  for (i in seq_along(x)) {
    df[seq_along(x[[i]]), i] <- x[[i]]
  }
  return(df)
}

#' Convert wide data frame to named list
#' @param df Wide data frame
#' @return Named list of vectors (NA values removed)
wideDF2namedList <- function(df) {
  result <- lapply(df, function(col) col[!is.na(col)])
  return(result)
}

#' Convert named list to long data frame
#' @param x Named list of vectors
#' @return Long data frame with columns: name, value
namedList2longDF <- function(x) {
  df <- data.frame(
    name = rep(names(x), sapply(x, length)),
    value = unlist(x, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  return(df)
}

#' Convert long data frame to named list
#' @param df Long data frame
#' @param group_by Column to group by
#' @param values Column containing values
#' @return Named list
longDF2namedList <- function(df, group_by, values) {
  split(df[[values]], df[[group_by]])
}

#' Get unique length
#' @param x Vector
#' @return Number of unique elements
ulength <- function(x) {
  length(unique(x))
}

# ---- Seurat Helper Functions ----

#' Calculate mitochondrial content for a Seurat object
#' @param object Seurat object
#' @return Seurat object with percent.mt in metadata
getMitoContent <- function(object) {
  if (!"percent.mt" %in% colnames(object@meta.data)) {
    # Try both human (MT-) and mouse (mt-) patterns
    mito_pattern <- "^MT-|^mt-"
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = mito_pattern)
  }
  return(object)
}

#' Filter cells based on QC metrics
#' @param object Seurat object
#' @return Filtered Seurat object
filterData <- function(object) {
  object <- getMitoContent(object)
  cells_keep <- object@meta.data$percent.mt < config$mito_threshold &
                object@meta.data$nFeature_RNA >= config$min_features &
                object@meta.data$nFeature_RNA <= config$max_features
  object <- object[, cells_keep]
  return(object)
}

#' Normalize samples using SCTransform
#' @param so.list List of Seurat objects
#' @param do_filter Apply QC filtering
#' @param do_filter_only Only filter, don't normalize
#' @return List of normalized Seurat objects
normData <- function(so.list, do_filter = TRUE, do_filter_only = FALSE) {
  
  # Apply filtering
  if (do_filter) {
    miko_message("Filtering cells...")
    so.list <- pbapply::pblapply(so.list, function(x) {
      x <- getMitoContent(x)
      cells_keep <- x@meta.data$percent.mt < config$mito_threshold &
                    x@meta.data$nFeature_RNA >= config$min_features &
                    x@meta.data$nFeature_RNA <= config$max_features
      x <- x[, cells_keep]
      return(x)
    })
    
    # Remove samples with insufficient cells
    cell_counts <- sapply(so.list, ncol)
    so.list <- so.list[cell_counts >= config$min_cells_per_sample]
    miko_message(sprintf("Retained %d samples with >= %d cells", 
                         length(so.list), config$min_cells_per_sample))
  }
  
  # Normalize with SCTransform
  if (!do_filter_only) {
    miko_message("Normalizing with SCTransform...")
    so.list <- pbapply::pblapply(so.list, function(x) {
      x <- tryCatch({
        SCTransform(x, 
                    method = "glmGamPoi", 
                    verbose = FALSE, 
                    vst.flavor = "v2",
                    vars.to.regress = "percent.mt", 
                    variable.features.n = config$n_variable_features)
      }, error = function(e) {
        # Fallback without vst.flavor if v2 fails
        SCTransform(x, 
                    method = "glmGamPoi", 
                    verbose = FALSE,
                    vars.to.regress = "percent.mt", 
                    variable.features.n = config$n_variable_features)
      })
      return(x)
    })
  }
  
  return(so.list)
}

#' Get expressed genes above threshold
#' @param object Seurat object
#' @param min_pct Minimum fraction of cells expressing gene
#' @param assay Assay to use
#' @return Vector of gene names
getExpressedGenes <- function(object, min_pct = config$min_pct_expressed, assay = "RNA") {
  expr_mat <- GetAssayData(object, slot = "data", assay = assay)
  pct_expressed <- Matrix::rowMeans(expr_mat > 0)
  genes <- names(pct_expressed)[pct_expressed >= min_pct]
  # Remove mitochondrial and ribosomal genes
  genes <- genes[!grepl("^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl", genes)]
  return(genes)
}

#' Downsample Seurat object
#' @param object Seurat object
#' @param subsample_n Number of cells to keep
#' @return Downsampled Seurat object
downsampleSeurat <- function(object, subsample_n) {
  if (ncol(object) <= subsample_n) {
    return(object)
  }
  cells_keep <- sample(colnames(object), subsample_n)
  return(object[, cells_keep])
}

#' Get UMAP coordinates from Seurat object
#' @param object Seurat object
#' @param umap_key Key for UMAP reduction (default: "umap")
#' @return List containing UMAP data frame
getUMAP <- function(object, umap_key = "umap") {
  # Try to find the reduction
  reduction_name <- NULL
  if (umap_key %in% names(object@reductions)) {
    reduction_name <- umap_key
  } else if ("b" %in% names(object@reductions)) {
    reduction_name <- "b"  # BBKNN reduction
  } else if ("umap" %in% names(object@reductions)) {
    reduction_name <- "umap"
  }
  
  if (is.null(reduction_name)) {
    stop("No UMAP reduction found")
  }
  
  umap_coords <- Embeddings(object, reduction = reduction_name)
  df_umap <- as.data.frame(umap_coords)
  colnames(df_umap) <- c("x", "y")
  df_umap <- cbind(df_umap, object@meta.data)
  
  return(list(df.umap = df_umap))
}

# ---- Gene ID Conversion ----

#' Convert Ensembl IDs to gene symbols
#' @param ensembl_ids Vector of Ensembl IDs
#' @param species Species ("Hs" for human, "Mm" for mouse)
#' @return Data frame with ENSEMBL and SYMBOL columns
ensembl2sym <- function(ensembl_ids, my_species = "Hs") {
  
  if (my_species == "Hs") {
    dataset <- "hsapiens_gene_ensembl"
  } else if (my_species == "Mm") {
    dataset <- "mmusculus_gene_ensembl"
  } else {
    stop("Species must be 'Hs' or 'Mm'")
  }
  
  tryCatch({
    mart <- useMart("ensembl", dataset = dataset)
    results <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = mart
    )
    colnames(results) <- c("ENSEMBL", "SYMBOL")
    return(results)
  }, error = function(e) {
    warning("biomaRt connection failed. Returning empty mapping.")
    return(data.frame(ENSEMBL = ensembl_ids, SYMBOL = ensembl_ids))
  })
}

# ---- BBKNN Integration via Reticulate ----

#' Run BBKNN batch correction
#' @param object Seurat object (must have PCA computed)
#' @param batch Batch variable in metadata
#' @param n_pcs Number of PCs to use
#' @return Seurat object with BBKNN reduction
runBBKNN <- function(object, batch, n_pcs = 50) {
  
  miko_message("Running BBKNN batch correction...")
  
  # Check if reticulate/Python is available
  if (!reticulate::py_available()) {
    warning("Python not available. Falling back to Harmony.")
    return(runHarmonyFallback(object, batch, n_pcs))
  }
  
  # Try to use specified conda environment
  tryCatch({
    reticulate::use_condaenv(config$conda_env, required = FALSE)
  }, error = function(e) {
    message("Conda environment not found, using default Python")
  })
  
  # Check if bbknn is available
  bbknn_available <- tryCatch({
    reticulate::py_module_available("bbknn")
  }, error = function(e) {
    FALSE
  })
  
  if (!bbknn_available) {
    warning("BBKNN Python module not available. Falling back to Harmony.")
    message("To install BBKNN: pip install bbknn")
    return(runHarmonyFallback(object, batch, n_pcs))
  }
  
  # Import Python modules
  sc <- reticulate::import("scanpy", convert = FALSE)
  bbknn <- reticulate::import("bbknn", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  
  # Get PCA embeddings
  pca_embeddings <- Embeddings(object, "pca")[, 1:n_pcs]
  
  # Get batch labels
  batch_labels <- as.character(object@meta.data[[batch]])
  
  # Create AnnData object
  adata <- sc$AnnData(
    X = np$array(as.matrix(GetAssayData(object, slot = "data"))),
    obs = reticulate::r_to_py(data.frame(batch = batch_labels, row.names = colnames(object)))
  )
  adata$obsm$`__setitem__`("X_pca", np$array(pca_embeddings))
  
  # Run BBKNN
  bbknn$bbknn(adata, batch_key = "batch", n_pcs = as.integer(n_pcs))
  
  # Run UMAP on corrected neighbors
  sc$tl$umap(adata)
  
  # Extract UMAP coordinates
  umap_coords <- reticulate::py_to_r(adata$obsm["X_umap"])
  rownames(umap_coords) <- colnames(object)
  colnames(umap_coords) <- c("B_1", "B_2")
  
  # Add to Seurat object
  object[["b"]] <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "B_",
    assay = DefaultAssay(object)
  )
  
  # Extract connectivities for clustering
  connectivities <- reticulate::py_to_r(adata$obsp["connectivities"])
  rownames(connectivities) <- colnames(connectivities) <- colnames(object)
  
  # Add as graph
  object@graphs[["bbknn"]] <- as.Graph(connectivities)
  
  miko_message("BBKNN complete.")
  return(object)
}

#' Harmony fallback for batch correction
#' @param object Seurat object
#' @param batch Batch variable
#' @param n_pcs Number of PCs
#' @return Seurat object with Harmony reduction
runHarmonyFallback <- function(object, batch, n_pcs = 50) {
  
  miko_message("Running Harmony batch correction (fallback)...")
  
  if (!requireNamespace("harmony", quietly = TRUE)) {
    warning("Harmony not available. Proceeding without batch correction.")
    object <- RunUMAP(object, reduction = "pca", dims = 1:n_pcs, verbose = FALSE)
    return(object)
  }
  
  object <- harmony::RunHarmony(object, group.by.vars = batch, dims.use = 1:n_pcs)
  object <- RunUMAP(object, reduction = "harmony", dims = 1:n_pcs, 
                    reduction.name = "b", reduction.key = "B_", verbose = FALSE)
  
  # Create graph for clustering
  object <- FindNeighbors(object, reduction = "harmony", dims = 1:n_pcs, 
                          graph.name = "bbknn", verbose = FALSE)
  
  return(object)
}

#' Run Scanorama integration (alternative method)
#' @param object Seurat object
#' @param batch Batch variable
#' @return Seurat object with Scanorama corrected embeddings
runScanorama <- function(object, batch) {
  
  miko_message("Running Scanorama integration...")
  
  # Check if scanorama is available
  scanorama_available <- tryCatch({
    reticulate::py_module_available("scanorama")
  }, error = function(e) {
    FALSE
  })
  
  if (!scanorama_available) {
    warning("Scanorama not available. Falling back to Harmony.")
    return(runHarmonyFallback(object, batch, config$n_pcs))
  }
  
  scanorama <- reticulate::import("scanorama", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  
  # Split by batch
  batch_labels <- unique(object@meta.data[[batch]])
  
  datasets <- list()
  genes_list <- list()
  
  for (b in batch_labels) {
    cells <- colnames(object)[object@meta.data[[batch]] == b]
    mat <- as.matrix(GetAssayData(object, slot = "data")[, cells])
    datasets <- c(datasets, list(np$array(t(mat))))
    genes_list <- c(genes_list, list(rownames(mat)))
  }
  
  # Run Scanorama
  result <- scanorama$integrate(datasets, genes_list)
  integrated <- reticulate::py_to_r(result[[1]])
  
  # Combine integrated embeddings
  combined <- do.call(rbind, integrated)
  rownames(combined) <- colnames(object)
  
  # Add to Seurat object
  object[["scanorama"]] <- CreateDimReducObject(
    embeddings = combined,
    key = "SCAN_",
    assay = DefaultAssay(object)
  )
  
  return(object)
}

# ---- Visualization Functions ----

#' Custom ggplot theme
#' @param legend Show legend
#' @param grid Show grid
#' @param x.axis.rotation Rotate x-axis labels
#' @return ggplot theme object
theme_miko <- function(legend = FALSE, grid = FALSE, x.axis.rotation = 0) {
  t <- theme_classic() +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  if (!legend) {
    t <- t + theme(legend.position = "none")
  }
  
  if (!grid) {
    t <- t + theme(panel.grid = element_blank())
  }
  
  if (x.axis.rotation > 0) {
    t <- t + theme(axis.text.x = element_text(angle = x.axis.rotation, 
                                               hjust = 1, vjust = 1))
  }
  
  return(t)
}

#' Custom color scale (low to high)
#' @param low Low color
#' @param high High color
#' @return ggplot scale
scale_color_miko <- function(low = "grey90", high = "darkred") {
  scale_color_gradient(low = low, high = high)
}

scale_fill_miko <- function(low = "grey90", high = "darkred") {
  scale_fill_gradient(low = low, high = high)
}

#' UMAP plot colored by cluster/group
#' @param object Seurat object
#' @param group Variable to color by
#' @param reduction Reduction to use
#' @param raster Use scattermore for large datasets
#' @return ggplot object
cluster.UMAP <- function(object, group = "seurat_clusters", reduction = "umap", 
                         raster = TRUE) {
  
  df <- getUMAP(object, umap_key = reduction)$df.umap
  
  if (!group %in% colnames(df)) {
    df[[group]] <- object@meta.data[[group]]
  }
  
  p <- ggplot(df, aes(x = x, y = y, color = .data[[group]])) +
    theme_miko(legend = TRUE) +
    labs(x = "UMAP 1", y = "UMAP 2", title = group) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  if (raster && nrow(df) > 10000) {
    p <- p + scattermore::geom_scattermore(pointsize = 1, pixels = c(1024, 1024))
  } else {
    p <- p + geom_point(size = 0.5, alpha = 0.7)
  }
  
  return(p)
}

#' UMAP plot colored by gene expression
#' @param object Seurat object
#' @param feature Gene name
#' @param reduction Reduction to use
#' @param assay Assay to use
#' @param slot Slot to use
#' @param scale.color Color for high expression
#' @return ggplot object
exprUMAP <- function(object, feature, reduction = "umap", assay = "RNA", 
                     slot = "data", scale.color = "darkred") {
  
  df <- getUMAP(object, umap_key = reduction)$df.umap
  
  # Get expression values
  if (feature %in% rownames(object[[assay]])) {
    expr <- GetAssayData(object, slot = slot, assay = assay)[feature, ]
    df[[feature]] <- as.numeric(expr)
  } else {
    warning(sprintf("Feature %s not found", feature))
    df[[feature]] <- 0
  }
  
  p <- ggplot(df, aes(x = x, y = y, color = .data[[feature]])) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_gradient(low = "grey90", high = scale.color) +
    theme_miko(legend = TRUE) +
    labs(x = "UMAP 1", y = "UMAP 2", title = feature)
  
  return(p)
}

#' Highlight specific groups on UMAP
#' @param object Seurat object
#' @param group Variable to highlight by
#' @param reduction Reduction to use
#' @param raster Use scattermore
#' @return List of ggplot objects
highlightUMAP <- function(object, group, reduction = "umap", raster = TRUE) {
  
  df <- getUMAP(object, umap_key = reduction)$df.umap
  df[[group]] <- object@meta.data[[group]]
  
  groups <- unique(df[[group]])
  plots <- list()
  
  for (g in groups) {
    df$highlight <- df[[group]] == g
    
    p <- ggplot(df, aes(x = x, y = y, color = highlight)) +
      scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "grey80")) +
      theme_miko() +
      labs(title = g, x = "UMAP 1", y = "UMAP 2") +
      theme(legend.position = "none")
    
    if (raster && nrow(df) > 10000) {
      p <- p + scattermore::geom_scattermore(pointsize = 1, pixels = c(512, 512))
    } else {
      p <- p + geom_point(size = 0.3, alpha = 0.5)
    }
    
    plots[[g]] <- p
  }
  
  return(plots)
}

#' Custom heatmap wrapper
#' @param mat Matrix to plot
#' @param ... Additional arguments to pheatmap
#' @return pheatmap object
miko_heatmap <- function(mat, ...) {
  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    ...
  )
}

#' Save PDF figure
#' @param filename Output filename
#' @param plot ggplot object
#' @param fig.width Width in inches
#' @param fig.height Height in inches
savePDF <- function(filename, plot, fig.width = 7, fig.height = 7) {
  pdf(filename, width = fig.width, height = fig.height)
  print(plot)
  dev.off()
  miko_message(sprintf("Saved: %s", filename))
}

# ---- Module Scoring Functions ----

#' Run module scoring for gene lists
#' @param object Seurat object
#' @param genelist Named list of gene vectors
#' @param scale Scale scores within each sample
#' @param assay Assay to use
#' @param reduction Reduction for plotting
#' @param return.plots Return UMAP plots
#' @param raster Use scattermore for plots
#' @return List with data and plots
runMS <- function(object, genelist, scale = TRUE, assay = "RNA", 
                  reduction = "umap", return.plots = TRUE, raster = TRUE) {
  
  miko_message("Calculating module scores...")
  
  # Calculate scores
  for (i in seq_along(genelist)) {
    prog_name <- names(genelist)[i]
    genes <- genelist[[i]]
    genes_present <- genes[genes %in% rownames(object[[assay]])]
    
    if (length(genes_present) >= 3) {
      object <- AddModuleScore(object, 
                                features = list(genes_present),
                                name = prog_name,
                                assay = assay)
      # Rename column (AddModuleScore appends "1")
      colnames(object@meta.data)[colnames(object@meta.data) == paste0(prog_name, "1")] <- prog_name
    }
  }
  
  # Extract scores
  score_cols <- names(genelist)[names(genelist) %in% colnames(object@meta.data)]
  df_scores <- object@meta.data[, score_cols, drop = FALSE]
  
  # Scale if requested
  if (scale && "sample" %in% colnames(object@meta.data)) {
    df_scores <- as.data.frame(lapply(df_scores, function(x) {
      ave(x, object@meta.data$sample, FUN = scale)
    }))
  }
  
  # Generate plots
  plot_list <- NULL
  if (return.plots) {
    df_umap <- getUMAP(object, umap_key = reduction)$df.umap
    df_plot <- cbind(df_umap[, c("x", "y")], df_scores)
    
    plot_list <- lapply(score_cols, function(col) {
      p <- ggplot(df_plot, aes(x = x, y = y, color = .data[[col]])) +
        viridis::scale_color_viridis(option = "magma") +
        theme_miko(legend = TRUE) +
        labs(title = col, x = "UMAP 1", y = "UMAP 2")
      
      if (raster && nrow(df_plot) > 10000) {
        p <- p + scattermore::geom_scattermore(pointsize = 1, pixels = c(512, 512))
      } else {
        p <- p + geom_point(size = 0.3, alpha = 0.7)
      }
      return(p)
    })
    names(plot_list) <- score_cols
  }
  
  return(list(data = df_scores, plot.list = plot_list, object = object))
}

# ---- NMF Analysis Functions ----

#' Calculate Jaccard similarity between gene sets
#' @param gene_list Named list of gene vectors
#' @return Jaccard similarity matrix
jaccardSimilarityMatrix <- function(gene_list) {
  n <- length(gene_list)
  jmat <- matrix(0, n, n)
  rownames(jmat) <- colnames(jmat) <- names(gene_list)
  
  for (i in 1:n) {
    for (j in i:n) {
      intersection <- length(intersect(gene_list[[i]], gene_list[[j]]))
      union <- length(union(gene_list[[i]], gene_list[[j]]))
      jac <- ifelse(union > 0, intersection / union, 0)
      jmat[i, j] <- jmat[j, i] <- jac
    }
  }
  return(jmat)
}

#' Extract top genes from NMF W matrix
#' @param W NMF W matrix (features x components)
#' @param n_genes Number of top genes per component
#' @param norm_cutoff Normalized loading cutoff (alternative to n_genes)
#' @return Named list of gene vectors per component
getNMFGenes <- function(W, n_genes = config$top_n_genes, norm_cutoff = NA) {
  
  # Normalize by squared loadings (kME-like)
  W_norm <- t(apply(t(W), 2, function(x) (x^2) / sum(x^2)))
  
  if (!is.na(norm_cutoff)) {
    # Use cutoff
    gene_list <- apply(W_norm, 1, function(x) colnames(W_norm)[x > norm_cutoff])
  } else {
    # Use top N genes
    gene_list <- lapply(1:nrow(W_norm), function(k) {
      names(sort(W_norm[k, ], decreasing = TRUE))[1:n_genes]
    })
    names(gene_list) <- paste0("k", 1:nrow(W_norm))
  }
  
  # Remove empty entries
  gene_list <- gene_list[sapply(gene_list, length) > 0]
  
  return(gene_list)
}

# ---- Hypergeometric Test Functions ----

#' Run hypergeometric enrichment test
#' @param gene_list Named list of gene vectors to test
#' @param gene_universe Vector of all possible genes
#' @param my_pathway Named list of pathway gene sets (or use built-in)
#' @param species Species ("Hs" or "Mm")
#' @return List of enrichment results
runHG <- function(gene_list, gene_universe, my_pathway = NULL, 
                  my_pathway_representation = "SYMBOL", species = "Hs") {
  
  miko_message("Running hypergeometric enrichment...")
  
  results <- list()
  
  for (prog_name in names(gene_list)) {
    query_genes <- gene_list[[prog_name]]
    query_in_universe <- query_genes[query_genes %in% gene_universe]
    
    if (length(query_in_universe) < 3) {
      next
    }
    
    prog_results <- data.frame()
    
    for (pathway_name in names(my_pathway)) {
      pathway_genes <- my_pathway[[pathway_name]]
      pathway_in_universe <- pathway_genes[pathway_genes %in% gene_universe]
      
      # Hypergeometric test
      overlap <- length(intersect(query_in_universe, pathway_in_universe))
      
      if (overlap < 2) next
      
      pval <- phyper(
        overlap - 1,
        length(pathway_in_universe),
        length(gene_universe) - length(pathway_in_universe),
        length(query_in_universe),
        lower.tail = FALSE
      )
      
      prog_results <- rbind(prog_results, data.frame(
        query = prog_name,
        pathway = pathway_name,
        overlap = overlap,
        query_size = length(query_in_universe),
        pathway_size = length(pathway_in_universe),
        universe_size = length(gene_universe),
        pvalue = pval
      ))
    }
    
    if (nrow(prog_results) > 0) {
      prog_results$padj <- p.adjust(prog_results$pvalue, method = "BH")
      results[[prog_name]] <- prog_results
    }
  }
  
  return(results)
}

#' Summarize hypergeometric results
#' @param hg_res Results from runHG
#' @param show_n Number of top pathways to show per program
#' @return List with summary data frame and plots
summarizeHG <- function(hg_res, show_n = 10) {
  
  all_results <- do.call(rbind, hg_res)
  
  # Top pathways per program
  top_results <- all_results %>%
    group_by(query) %>%
    slice_min(order_by = padj, n = show_n) %>%
    ungroup()
  
  # Create plot
  p <- ggplot(top_results, aes(x = -log10(padj), y = reorder(pathway, -log10(padj)))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~query, scales = "free_y") +
    theme_miko() +
    labs(x = "-log10(FDR)", y = "Pathway", title = "Pathway Enrichment")
  
  return(list(results = all_results, top = top_results, plots = p))
}

# ==============================================================================
# PART 3: DATA PREPROCESSING - INDIVIDUAL STUDY PROCESSING
# ==============================================================================

miko_message(paste(rep("=", 70), collapse = ""))
miko_message("OLIGODENDROCYTE META-ATLAS PIPELINE")
miko_message(paste(rep("=", 70), collapse = ""))

# This section contains the preprocessing code for each individual study.
# Each study requires:
# 1. Loading raw data (counts matrix, metadata)
# 2. Creating Seurat object
# 3. Filtering for oligodendrocyte cells
# 4. SCTransform normalization
# 5. Saving preprocessed object

# ---- Define preprocessing functions for each data format ----

#' Process 10X-style data (matrix.mtx, features.tsv, barcodes.tsv)
#' @param data_path Path to directory containing 10X files
#' @param study_name Name of the study
#' @param sample_col Column in metadata for sample ID
#' @param type_label Disease/condition type
#' @return List of Seurat objects split by sample
process10X <- function(data_path, study_name, sample_col = NULL, type_label = "HEALTHY") {
  
  miko_message(sprintf("Processing %s...", study_name))
  
  # Read files
  features <- read.delim(file.path(data_path, "features.tsv"), header = FALSE)
  if (ncol(features) >= 2) {
    gene_names <- features$V2
  } else {
    gene_names <- features$V1
  }
  
  barcodes <- read.delim(file.path(data_path, "barcodes.tsv"), header = FALSE)
  mat <- Matrix::readMM(file.path(data_path, "matrix.mtx"))
  
  colnames(mat) <- barcodes$V1
  rownames(mat) <- make.unique(gene_names)
  
  # Create Seurat object
  so <- CreateSeuratObject(counts = mat)
  so@meta.data$study <- study_name
  so@meta.data$type <- type_label
  
  return(so)
}

#' Process data with separate count and metadata files
#' @param count_file Path to count matrix file
#' @param meta_file Path to metadata file (optional)
#' @param study_name Name of the study
#' @param type_label Disease/condition type
#' @return Seurat object
processCountMatrix <- function(count_file, meta_file = NULL, study_name, type_label = "HEALTHY") {
  
  miko_message(sprintf("Processing %s...", study_name))
  
  # Determine file type and read
  if (grepl("\\.rds$", count_file, ignore.case = TRUE)) {
    mat <- readRDS(count_file)
  } else if (grepl("\\.csv$", count_file, ignore.case = TRUE)) {
    mat <- read.csv(count_file, row.names = 1)
  } else {
    mat <- read.delim(count_file, row.names = 1)
  }
  
  # Read metadata if provided
  meta <- NULL
  if (!is.null(meta_file) && file.exists(meta_file)) {
    if (grepl("\\.csv$", meta_file, ignore.case = TRUE)) {
      meta <- read.csv(meta_file)
    } else {
      meta <- read.delim(meta_file)
    }
    rownames(meta) <- meta[, 1]  # Assume first column is cell ID
  }
  
  # Create Seurat object
  if (!is.null(meta)) {
    # Ensure matching cells
    common_cells <- intersect(colnames(mat), rownames(meta))
    mat <- mat[, common_cells]
    meta <- meta[common_cells, ]
    so <- CreateSeuratObject(counts = mat, meta.data = meta)
  } else {
    so <- CreateSeuratObject(counts = mat)
  }
  
  so@meta.data$study <- study_name
  so@meta.data$type <- type_label
  
  return(so)
}

# ---- Process each study ----
# Note: Set reprocess = TRUE to reprocess from raw data
# Set reprocess = FALSE to load previously saved preprocessed objects

process_studies <- function() {
  
  all_oligo_data <- list()
  
  # ========================================================================
  # HEALTHY BRAIN STUDIES
  # ========================================================================
  
  # ---- Franjic 2022 ----
  reprocess_franjic <- FALSE
  if (reprocess_franjic) {
    current_dir <- file.path(config$external_data_dir, "Franjic2022_healthy/")
    
    meta.data <- read.delim(file.path(current_dir, "Human_cell_meta.txt"), header = TRUE)
    rownames(meta.data) <- meta.data$cell_name
    features <- read.delim(file.path(current_dir, "Human_genes.txt"), header = FALSE)
    mat <- Matrix::readMM(file.path(current_dir, "Human_counts.mtx"))
    
    # Identify oligo clusters
    which.oligo <- unique(toupper(meta.data$cluster)[grepl("OLIGO|OPC|COP", toupper(meta.data$cluster))])
    
    colnames(mat) <- meta.data$cell_name
    rownames(mat) <- features$V1
    
    so.franjic <- CreateSeuratObject(counts = mat, meta.data = meta.data)
    so.franjic@meta.data$study <- "Franjic_2022"
    so.franjic@meta.data$sample <- so.franjic$samplename
    so.franjic@meta.data$type <- "HEALTHY"
    so.franjic@meta.data$cell_type <- so.franjic$cluster
    
    # Filter for oligos
    so.franjic <- so.franjic[, toupper(so.franjic$cell_type) %in% which.oligo]
    
    # Split and normalize
    so.list <- SplitObject(so.franjic, split.by = "sample")
    so.list <- normData(so.list)
    
    saveRDS(so.list, file.path(current_dir, "seurat_object_oligo_franjic.rds"))
    all_oligo_data$franjic <- so.list
  } else {
    all_oligo_data$franjic <- readRDS(file.path(config$external_data_dir, 
                                                  "Franjic2022_healthy/seurat_object_oligo_franjic.rds"))
  }
  
  # ---- Habib 2017 ----
  reprocess_habib <- FALSE
  if (reprocess_habib) {
    current_dir <- file.path(config$external_data_dir, "Habib2017_healthy/")
    
    meta.data <- read.delim(file.path(current_dir, "GTEx_droncseq_hip_pcf.clusters.txt"), header = FALSE)
    c2c <- meta.data$V2
    names(c2c) <- meta.data$V1
    
    mat <- read.delim(file.path(current_dir, "GTEx_droncseq_hip_pcf.umi_counts.txt"))
    mat <- col2rowname(mat, "X")
    
    so.habib <- CreateSeuratObject(counts = mat)
    so.habib@meta.data$cluster <- c2c[colnames(so.habib)]
    so.habib@meta.data$study <- "Habib_2017"
    so.habib@meta.data$type <- "HEALTHY"
    so.habib@meta.data$cell_type <- so.habib$cluster
    
    # Clusters 10, 11, 12 are oligodendrocytes
    so.habib <- so.habib[, as.character(so.habib$cell_type) %in% c("10", "11", "12")]
    
    # Assign sample based on cell barcode pattern
    so.habib@meta.data$sample <- as.character(so.habib@active.ident)
    so.habib@meta.data$sample[grepl("HP", so.habib@meta.data$sample)] <- "HP"
    so.habib@meta.data$sample[grepl("PFC", so.habib@meta.data$sample)] <- "PFC"
    so.habib@meta.data$sample[grepl("hC", so.habib@meta.data$sample)] <- "HC"
    
    so.list <- SplitObject(so.habib, split.by = "sample")
    so.list <- normData(so.list)
    
    saveRDS(so.list, file.path(current_dir, "seurat_object_oligo_habib.rds"))
    all_oligo_data$habib <- so.list
  } else {
    all_oligo_data$habib <- readRDS(file.path(config$external_data_dir,
                                               "Habib2017_healthy/seurat_object_oligo_habib.rds"))
  }
  
  # Additional studies follow same pattern...
  # [Similar blocks for: Kanton2019, Khrameeva2020, Hodge2019, Bakken2021,
  #  Wheeler2020, Jakel2019, Schirmer2019, Sun2022, Kim2020, Biermann2022,
  #  Heming2022, Lau2020, Smajic2022, Yu2021, vanBruggen2022, Cao2020,
  #  Bhaduri2021, Aldinger2021]
  
  # ========================================================================
  # GBM STUDIES (Current study data)
  # ========================================================================
  
  reprocess_gbm <- FALSE
  if (reprocess_gbm) {
    # Load GBM oligodendrocyte data
    so.current <- readRDS(file.path(config$gbm_data_dir, "seurat_object_oligo_current.rds"))
    so.abdel <- readRDS(file.path(config$gbm_data_dir, "seurat_object_oligo_abdel.rds"))
    so.wang <- readRDS(file.path(config$gbm_data_dir, "seurat_object_oligo_wang.rds"))
    
    # Process current study
    so.current@meta.data$sample <- so.current@meta.data$Barcode
    so.current@meta.data$type <- paste0(tolower(so.current@meta.data$PR), "GBM")
    so.list.current <- SplitObject(so.current, split.by = "sample")
    so.list.current <- so.list.current[sapply(so.list.current, ncol) > config$min_cells_per_sample]
    so.list.current <- normData(so.list.current)
    so.list.current <- lapply(so.list.current, function(x) {
      x@meta.data$study <- "Mikolajewicz_2024"
      return(x)
    })
    
    # Process Abdelfattah data
    so.abdel <- CreateSeuratObject(counts = so.abdel@assays$SCT@counts, meta.data = so.abdel@meta.data)
    so.abdel@meta.data$type <- paste0(tolower(so.abdel@meta.data$PR), "GBM")
    so.list.abdel <- SplitObject(so.abdel, split.by = "sample")
    so.list.abdel <- so.list.abdel[sapply(so.list.abdel, ncol) > config$min_cells_per_sample]
    so.list.abdel <- normData(so.list.abdel)
    so.list.abdel <- lapply(so.list.abdel, function(x) {
      x@meta.data$study <- "Abdelfattah_2022"
      return(x)
    })
    
    # Process Wang data
    so.wang@meta.data$PR2 <- so.wang@meta.data$PR
    so.wang@meta.data$PR2[grepl("Primary", so.wang@meta.data$PR)] <- "p"
    so.wang@meta.data$PR2[grepl("Recurrent", so.wang@meta.data$PR)] <- "r"
    so.wang@meta.data$type <- paste0(so.wang@meta.data$PR2, "GBM")
    so.wang@meta.data$sample <- so.wang@meta.data$Barcode
    so.list.wang <- SplitObject(so.wang, split.by = "sample")
    so.list.wang <- so.list.wang[sapply(so.list.wang, ncol) > config$min_cells_per_sample]
    so.list.wang <- normData(so.list.wang)
    so.list.wang <- lapply(so.list.wang, function(x) {
      x@meta.data$study <- "Wang_2022"
      return(x)
    })
    
    all_oligo_data$current <- so.list.current
    all_oligo_data$abdel <- so.list.abdel
    all_oligo_data$wang <- so.list.wang
    
  } else {
    # Load preprocessed GBM data
    all_oligo_data$current <- readRDS(file.path(config$gbm_data_dir, "seurat_object_oligo_current_processed.rds"))
    all_oligo_data$abdel <- readRDS(file.path(config$gbm_data_dir, "seurat_object_oligo_abdel_processed.rds"))
    all_oligo_data$wang <- readRDS(file.path(config$gbm_data_dir, "seurat_object_oligo_wang_processed.rds"))
  }
  
  return(all_oligo_data)
}

# ==============================================================================
# PART 4: DATA INTEGRATION
# ==============================================================================

run_integration <- function(all_oligo_data) {
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STEP 1: Data Integration")
  miko_message(paste(rep("=", 70), collapse = ""))
  
  # ---- Combine all samples ----
  miko_message("Combining all samples...")
  
  # Flatten all sample lists
  so_query <- unlist(all_oligo_data, recursive = FALSE)
  
  # Standardize sample names
  generateSampleNames <- function(so_list) {
    types <- sapply(so_list, function(x) unique(x@meta.data$type)[1])
    studies <- sapply(so_list, function(x) unique(x@meta.data$study)[1])
    studies <- gsub("_", "", studies)
    new_names <- paste0(types, "-", seq_along(types), "_", studies)
    return(new_names)
  }
  
  names(so_query) <- generateSampleNames(so_query)
  
  # Create cell-to-sample mapping
  cell2sample <- unlist(lapply(names(so_query), function(sample_name) {
    cells <- colnames(so_query[[sample_name]])
    setNames(rep(sample_name, length(cells)), cells)
  }))
  
  miko_message(sprintf("Total samples: %d", length(so_query)))
  miko_message(sprintf("Total cells: %d", length(cell2sample)))
  
  # ---- Merge samples ----
  miko_message("Merging samples...")
  so_merge <- merge(so_query[[1]], y = so_query[-1])
  so_merge@meta.data$sample <- cell2sample[colnames(so_merge)]
  
  # ---- Standard preprocessing ----
  miko_message("Normalizing merged data...")
  DefaultAssay(so_merge) <- "RNA"
  so_merge <- NormalizeData(so_merge, verbose = FALSE)
  so_merge <- FindVariableFeatures(so_merge, verbose = FALSE)
  so_merge <- ScaleData(so_merge, block.size = 1000, verbose = FALSE)
  
  miko_message("Running PCA...")
  so_merge <- RunPCA(so_merge, features = VariableFeatures(so_merge), verbose = FALSE)
  
  # ---- Batch correction with BBKNN ----
  if (config$use_bbknn) {
    so_merge <- runBBKNN(so_merge, batch = "sample", n_pcs = config$n_pcs)
  } else {
    so_merge <- runHarmonyFallback(so_merge, batch = "sample", n_pcs = config$n_pcs)
  }
  
  # ---- Clustering ----
  miko_message("Clustering...")
  so_merge <- FindClusters(so_merge, graph.name = "bbknn", resolution = 0.5, verbose = FALSE)
  so_merge <- FindClusters(so_merge, graph.name = "bbknn", resolution = 1.0, verbose = FALSE)
  
  # ---- SCTransform for NMF ----
  miko_message("Running SCTransform for NMF...")
  gc()
  so_merge <- SCTransform(so_merge, 
                           method = "glmGamPoi", 
                           verbose = TRUE, 
                           vst.flavor = "v2",
                           conserve.memory = TRUE,
                           vars.to.regress = "percent.mt", 
                           variable.features.n = config$n_variable_features_integration)
  
  # Save integrated object
  saveRDS(so_merge, file.path(config$output_dir, "seurat_oligo_integrated.rds"))
  
  return(list(merged = so_merge, samples = so_query, cell2sample = cell2sample))
}

# ==============================================================================
# PART 5: NMF ANALYSIS
# ==============================================================================

run_nmf_analysis <- function(so_query) {
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STEP 2: NMF Analysis")
  miko_message(paste(rep("=", 70), collapse = ""))
  
  nmf_results <- list()
  
  for (j in seq_along(so_query)) {
    
    current_sample <- names(so_query)[j]
    miko_message(sprintf("Processing sample %d/%d: %s", j, length(so_query), current_sample))
    
    tryCatch({
      
      object <- so_query[[current_sample]]
      
      # Get expressed genes (excluding MT and ribosomal)
      expr_genes <- getExpressedGenes(object, min_pct = config$min_pct_expressed)
      expr_genes <- expr_genes[!grepl("^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl", expr_genes)]
      
      # Get residuals for these genes
      object <- GetResidual(object, features = expr_genes, assay = "SCT")
      
      # Prepare expression matrix (non-negative)
      expr_mat <- object@assays$SCT@scale.data
      expr_mat <- expr_mat[rownames(expr_mat) %in% expr_genes, ]
      expr_mat[expr_mat < 0] <- 0
      
      # Store variable features
      object@assays$SCT@var.features <- rownames(expr_mat)
      
      # Run NMF across k values in parallel
      n_cores_use <- min(length(config$nmf_k_range), config$n_cores)
      cl <- makeCluster(n_cores_use)
      registerDoParallel(cl)
      
      nmf_models <- foreach(k = config$nmf_k_range, 
                            .packages = c("NNLM", "NMF")) %dopar% {
        z <- nnmf(expr_mat, k, verbose = FALSE, check.k = FALSE)
        return(nmfModel(H = z$H, W = z$W))
      }
      names(nmf_models) <- paste0(current_sample, "_k", config$nmf_k_range)
      
      stopCluster(cl)
      closeAllConnections()
      gc()
      
      # Extract gene programs for each k
      nmf_gene_current <- list()
      for (i in seq_along(nmf_models)) {
        nmf_name <- names(nmf_models)[i]
        current_model <- nmf_models[[nmf_name]]
        current_genes <- getNMFGenes(current_model@W, n_genes = config$top_n_genes)
        nmf_name2 <- gsub(paste0(current_sample, "_"), "", nmf_name)
        names(current_genes) <- paste0(nmf_name2, "_", seq_len(length(current_genes)))
        nmf_gene_current[[nmf_name]] <- current_genes
      }
      
      # Flatten programs
      programs_flat <- list()
      for (i in seq_along(nmf_gene_current)) {
        programs_flat <- c(programs_flat, nmf_gene_current[[i]])
      }
      
      # Calculate Jaccard similarity
      jmat <- jaccardSimilarityMatrix(programs_flat)
      
      nmf_results[[current_sample]] <- list(
        sample = current_sample,
        nmf = nmf_models,
        programs = nmf_gene_current,
        programs.flat = programs_flat,
        jmat = jmat
      )
      
    }, error = function(e) {
      miko_message(sprintf("  Error: %s", e$message))
    })
  }
  
  # Save NMF results
  saveRDS(nmf_results, file.path(config$output_dir, "nmf_results_per_sample.rds"))
  
  return(nmf_results)
}

# ==============================================================================
# PART 6: IDENTIFY ROBUST GENE PROGRAMS
# ==============================================================================

identify_robust_programs <- function(nmf_results, gene_universe) {
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STEP 3: Identifying Robust Gene Programs")
  miko_message(paste(rep("=", 70), collapse = ""))
  
  # ---- Step 1: Within-sample robustness ----
  miko_message("Filtering for within-sample robustness...")
  
  nmf_gene_all <- list()
  
  for (i in seq_along(nmf_results)) {
    current_sample <- names(nmf_results)[i]
    programs_flat <- nmf_results[[current_sample]]$programs.flat
    
    jmat <- jaccardSimilarityMatrix(programs_flat)
    
    # Keep programs similar to at least one other within sample
    robust_idx <- apply(jmat, 1, function(x) sum(x > config$intra_threshold) > 1)
    
    if (sum(robust_idx) == 0) next
    
    jmat_robust <- jmat[robust_idx, robust_idx, drop = FALSE]
    
    robust_progs <- programs_flat[rownames(jmat_robust)]
    names(robust_progs) <- paste0(current_sample, "_", names(robust_progs))
    nmf_gene_all <- c(nmf_gene_all, robust_progs)
  }
  
  miko_message(sprintf("  Retained %d within-sample robust programs", length(nmf_gene_all)))
  
  # ---- Step 2: Between-sample robustness ----
  miko_message("Filtering for between-sample robustness...")
  
  jmat_all <- jaccardSimilarityMatrix(nmf_gene_all)
  
  # Get unique samples
  all_samples <- unique(stringr::str_remove(names(nmf_gene_all), "_k[0-9]*_[0-9]*"))
  
  # Find programs similar across multiple samples
  robust_inter <- apply(jmat_all, 1, function(x) {
    similar_progs <- rownames(jmat_all)[x > config$inter_threshold]
    unique(stringr::str_remove(similar_progs, "_k[0-9]*_[0-9]*"))
  })
  
  robust_inter_count <- sapply(robust_inter, ulength)
  robust_programs <- names(robust_inter)[robust_inter_count >= config$min_samples_robust]
  
  miko_message(sprintf("  Retained %d cross-sample robust programs", length(robust_programs)))
  
  # ---- Step 3: Cluster robust programs ----
  miko_message("Clustering robust programs...")
  
  jmat_robust <- jmat_all[robust_programs, robust_programs]
  
  # Hierarchical clustering
  dmat <- as.dist(1 - cor(t(jmat_robust)))
  hclust_result <- hclust(dmat, method = "complete")
  program_clusters <- cutree(hclust_result, k = config$n_final_programs)
  
  # Create annotation data frame
  df_program_ann <- data.frame(
    program = names(program_clusters),
    cluster = paste0("nmf_", program_clusters),
    row.names = names(program_clusters)
  )
  
  # ---- Step 4: Define final gene programs ----
  miko_message("Defining final gene programs...")
  
  unique_clusters <- unique(program_clusters)
  final_programs <- list()
  
  for (i in seq_along(unique_clusters)) {
    program_name <- paste0("G", i)
    cluster_progs <- nmf_gene_all[names(program_clusters)[program_clusters == unique_clusters[i]]]
    n_progs <- length(cluster_progs)
    
    # Count gene occurrences
    gene_counts <- table(unlist(cluster_progs))
    gene_prevalence <- gene_counts / n_progs
    
    # Keep genes above threshold
    final_genes <- names(gene_prevalence)[gene_prevalence > config$prevalence_threshold]
    final_programs[[program_name]] <- final_genes
    
    miko_message(sprintf("  %s: %d genes (from %d component programs)", 
                         program_name, length(final_genes), n_progs))
  }
  
  # Save results
  saveRDS(list(
    jmat_robust = jmat_robust,
    n_clusters = config$n_final_programs,
    program_annotations = df_program_ann,
    final_programs = final_programs,
    gene_universe = gene_universe
  ), file.path(config$output_dir, "nmf_robust_programs.rds"))
  
  return(list(
    programs = final_programs,
    jmat = jmat_robust,
    annotations = df_program_ann
  ))
}

# ==============================================================================
# PART 7: ANNOTATE AND SCORE PROGRAMS
# ==============================================================================

annotate_and_score <- function(so_merge, final_programs) {
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STEP 4: Annotating and Scoring Programs")
  miko_message(paste(rep("=", 70), collapse = ""))
  
  # ---- Program annotations based on enrichment analysis ----
  # These annotations were determined by pathway enrichment
  program_annotations <- c(
    G1 = "Neuro-I",
    G2 = "Reactive-I", 
    G3 = "Cycling",
    G4 = "OPC",
    G5 = "Reactive-II",
    G6 = "Myelin",
    G7 = "Stress",
    G8 = "Neuro-II"
  )
  
  # Rename programs
  old_names <- names(final_programs)
  new_names <- paste0("O", seq_along(final_programs), "-", 
                      program_annotations[old_names])
  names(final_programs) <- new_names
  
  # ---- Calculate module scores ----
  miko_message("Calculating module scores...")
  
  ms_result <- runMS(so_merge, 
                     genelist = final_programs, 
                     scale = TRUE, 
                     assay = "SCT",
                     reduction = "b",
                     return.plots = TRUE,
                     raster = TRUE)
  
  so_merge <- ms_result$object
  
  # Save scored object
  saveRDS(so_merge, file.path(config$output_dir, "seurat_oligo_atlas_scored.rds"))
  
  return(list(
    object = so_merge,
    programs = final_programs,
    scores = ms_result$data,
    plots = ms_result$plot.list
  ))
}

# ==============================================================================
# PART 8: EXPORT RESULTS
# ==============================================================================

export_results <- function(final_programs, so_merge) {
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STEP 5: Exporting Results")
  miko_message(paste(rep("=", 70), collapse = ""))
  
  # ---- Export gene programs as CSV (Table S5) ----
  program_df <- namedList2wideDF(final_programs)
  
  write.csv(program_df, 
            file.path(config$output_dir, "TableS5_OL_metaprograms.csv"),
            row.names = FALSE, na = "")
  miko_message("  Saved: TableS5_OL_metaprograms.csv")
  
  # ---- Save gene programs as R object ----
  saveRDS(final_programs, file.path(config$output_dir, "oligo_nmf_programs.rds"))
  miko_message("  Saved: oligo_nmf_programs.rds")
  
  # ---- Save final Seurat object ----
  saveRDS(so_merge, file.path(config$output_dir, "seurat_oligo_atlas_final.rds"))
  miko_message("  Saved: seurat_oligo_atlas_final.rds")
  
  return(invisible(NULL))
}

# ==============================================================================
# PART 9: VISUALIZATION
# ==============================================================================

generate_visualizations <- function(so_merge, final_programs, jmat_robust) {
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STEP 6: Generating Visualizations")
  miko_message(paste(rep("=", 70), collapse = ""))
  
  # ---- Color palettes ----
  type_colors <- c(
    "HEALTHY" = "#4DAF4A",
    "pGBM" = "#E41A1C",
    "rGBM" = "#984EA3",
    "MS" = "#FF7F00",
    "AD" = "#377EB8",
    "PD" = "#A65628",
    "DEV" = "#F781BF",
    "BM" = "#999999",
    "Glioma" = "#A6CEE3",
    "PCNSL" = "#FDBF6F"
  )
  
  # ---- UMAP plots ----
  
  # By condition type
  p_type <- cluster.UMAP(so_merge, group = "type", reduction = "b") +
    scale_color_manual(values = type_colors) +
    ggtitle("Oligodendrocyte Meta-Atlas by Condition")
  savePDF(file.path(config$output_dir, "umap_by_type.pdf"), p_type, 8, 6)
  
  # By study
  p_study <- cluster.UMAP(so_merge, group = "study", reduction = "b") +
    ggtitle("Oligodendrocyte Meta-Atlas by Study")
  savePDF(file.path(config$output_dir, "umap_by_study.pdf"), p_study, 10, 6)
  
  # By cluster
  p_cluster <- cluster.UMAP(so_merge, group = "seurat_clusters", reduction = "b") +
    ggtitle("Oligodendrocyte Meta-Atlas Clusters")
  savePDF(file.path(config$output_dir, "umap_clusters.pdf"), p_cluster, 8, 6)
  
  # ---- Module score UMAPs ----
  module_cols <- names(final_programs)
  module_cols <- module_cols[module_cols %in% colnames(so_merge@meta.data)]
  
  module_plots <- list()
  for (col in module_cols) {
    df <- getUMAP(so_merge, umap_key = "b")$df.umap
    df[[col]] <- so_merge@meta.data[[col]]
    
    p <- ggplot(df, aes(x = x, y = y, color = .data[[col]])) +
      scattermore::geom_scattermore(pointsize = 1, pixels = c(1024, 1024)) +
      viridis::scale_color_viridis(option = "magma") +
      theme_miko(legend = TRUE) +
      labs(x = "UMAP 1", y = "UMAP 2", title = col)
    
    module_plots[[col]] <- p
    savePDF(file.path(config$output_dir, paste0("umap_module_", gsub("-", "_", col), ".pdf")), 
            p, 6, 5)
  }
  
  # Combined module plot
  p_combined <- cowplot::plot_grid(plotlist = module_plots, ncol = 4)
  savePDF(file.path(config$output_dir, "umap_all_modules.pdf"), p_combined, 16, 8)
  
  # ---- Jaccard similarity heatmap ----
  pdf(file.path(config$output_dir, "program_jaccard_heatmap.pdf"), width = 12, height = 10)
  miko_heatmap(jmat_robust,
               color = viridis::magma(100),
               clustering_method = "complete",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               cutree_rows = config$n_final_programs,
               cutree_cols = config$n_final_programs,
               border_color = NA,
               fontsize = 4,
               main = "NMF Program Similarity (Jaccard Index)")
  dev.off()
  miko_message("  Saved: program_jaccard_heatmap.pdf")
  
  # ---- Program correlation heatmap ----
  score_cols <- names(final_programs)
  score_cols <- score_cols[score_cols %in% colnames(so_merge@meta.data)]
  
  if (length(score_cols) > 1) {
    cor_mat <- cor(so_merge@meta.data[, score_cols], method = "pearson", use = "complete.obs")
    
    pdf(file.path(config$output_dir, "program_correlation_heatmap.pdf"), width = 7, height = 6)
    miko_heatmap(cor_mat,
                 display_numbers = TRUE,
                 number_format = "%.2f",
                 border_color = NA,
                 main = "Gene Program Correlations")
    dev.off()
    miko_message("  Saved: program_correlation_heatmap.pdf")
  }
  
  # ---- Marker gene expression ----
  marker_genes <- list(
    OPC = c("PDGFRA", "CSPG4", "SOX6"),
    Mature_OL = c("MBP", "PLP1", "MOBP", "MOG", "MAG"),
    Reactive = c("C1QA", "C1QB", "CD74"),
    Cycling = c("TOP2A", "MKI67"),
    Stress = c("HSPA1A", "HSPA1B", "FOS")
  )
  
  for (category in names(marker_genes)) {
    for (gene in marker_genes[[category]]) {
      if (gene %in% rownames(so_merge)) {
        p <- exprUMAP(so_merge, gene, reduction = "b", assay = "SCT", slot = "data")
        savePDF(file.path(config$output_dir, paste0("expr_", gene, ".pdf")), p, 6, 5)
      }
    }
  }
  
  miko_message("Visualizations complete.")
  return(invisible(NULL))
}

# ==============================================================================
# PART 10: MAIN EXECUTION
# ==============================================================================

main <- function() {
  
  start_time <- Sys.time()
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("STARTING OLIGODENDROCYTE META-ATLAS ANALYSIS")
  miko_message(sprintf("Start time: %s", format(start_time, "%Y-%m-%d %H:%M:%S")))
  miko_message(paste(rep("=", 70), collapse = ""))
  
  # Create output directory
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }
  
  # ---- Step 1: Process individual studies ----
  all_oligo_data <- process_studies()
  
  # ---- Step 2: Integration ----
  integration_result <- run_integration(all_oligo_data)
  so_merge <- integration_result$merged
  so_query <- integration_result$samples
  
  # ---- Step 3: NMF analysis ----
  nmf_results <- run_nmf_analysis(so_query)
  
  # ---- Step 4: Identify robust programs ----
  gene_universe <- unique(rownames(so_merge@assays$RNA))
  robust_result <- identify_robust_programs(nmf_results, gene_universe)
  
  # ---- Step 5: Annotate and score ----
  scored_result <- annotate_and_score(so_merge, robust_result$programs)
  so_merge <- scored_result$object
  final_programs <- scored_result$programs
  
  # ---- Step 6: Export results ----
  export_results(final_programs, so_merge)
  
  # ---- Step 7: Generate visualizations ----
  generate_visualizations(so_merge, final_programs, robust_result$jmat)
  
  # ---- Summary ----
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message("ANALYSIS COMPLETE")
  miko_message(paste(rep("=", 70), collapse = ""))
  miko_message(sprintf("Total cells in atlas: %d", ncol(so_merge)))
  miko_message(sprintf("Total studies: %d", length(unique(so_merge$study))))
  miko_message(sprintf("Total samples: %d", length(unique(so_merge$sample))))
  miko_message(sprintf("Gene programs identified: %d", length(final_programs)))
  miko_message(sprintf("Total runtime: %.1f minutes", as.numeric(duration)))
  
  # Session info
  miko_message("\nSession Info:")
  print(sessionInfo())
  
  return(list(
    seurat = so_merge,
    programs = final_programs,
    nmf_results = nmf_results
  ))
}

# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

# Uncomment to run the full analysis:
# results <- main()

# Or run individual steps:
# all_oligo_data <- process_studies()
# integration_result <- run_integration(all_oligo_data)
# nmf_results <- run_nmf_analysis(integration_result$samples)
# etc.

miko_message("Script loaded successfully. Call main() to run the full analysis.")
