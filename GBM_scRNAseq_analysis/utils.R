# ==============================================================================
# utils.R - Core Utility Functions
# ==============================================================================
# General utility functions for the OL-GBM analysis pipeline.
# ==============================================================================

#' Load configuration from YAML file
#'
#' @param config_path Path to config.yml file
#' @return Named list with configuration parameters
#' @export
load_config <- function(config_path = "config/config.yml") {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found: %s", config_path))
  }
  
  config <- yaml::read_yaml(config_path)
  
  # Set seed for reproducibility
  set.seed(config$reproducibility$seed)
  
  message(sprintf("[%s] Configuration loaded from: %s", 
                  Sys.time(), config_path))
  
  return(config)
}

#' Create output directories if they don't exist
#'
#' @param config Configuration list with paths
#' @export
setup_directories <- function(config) {
  dirs_to_create <- c(
    config$paths$results_dir,
    config$paths$figures_dir,
    config$paths$tables_dir,
    config$paths$processed_data_dir
  )
  
  for (dir_path in dirs_to_create) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("Created directory: %s", dir_path))
    }
  }
}

#' Log message with timestamp
#'
#' @param msg Message to log
#' @param level Log level: "INFO", "WARNING", "ERROR"
#' @export
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

#' Validate Seurat object has required metadata columns
#'
#' @param object Seurat object
#' @param required_cols Character vector of required column names
#' @return TRUE if valid, stops with error otherwise
#' @export
validate_seurat_metadata <- function(object, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(object@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required metadata columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  
  return(TRUE)
}

#' Check if required genes are present in Seurat object
#'
#' @param object Seurat object
#' @param genes Character vector of gene names
#' @param min_fraction Minimum fraction of genes that must be present (0-1)
#' @return Genes present in object
#' @export
check_genes_present <- function(object, genes, min_fraction = 0.5) {
  all_genes <- rownames(object)
  present_genes <- intersect(genes, all_genes)
  
  fraction_present <- length(present_genes) / length(genes)
  
  if (fraction_present < min_fraction) {
    warning(sprintf("Only %.1f%% of requested genes found in object",
                    fraction_present * 100))
  }
  
  log_message(sprintf("Found %d/%d genes (%.1f%%)", 
                      length(present_genes), 
                      length(genes),
                      fraction_present * 100))
  
  return(present_genes)
}

#' Convert long data.frame to named list
#'
#' @param df Data frame with group and value columns
#' @param group_by Column name for grouping
#' @param values Column name for values
#' @return Named list
#' @export
longDF_to_named_list <- function(df, group_by, values) {
  split(df[[values]], df[[group_by]])
}

#' Convert wide data.frame to named list (columns become list elements)
#'
#' @param df Wide data frame
#' @return Named list where each column is an element
#' @export
wideDF_to_named_list <- function(df) {
  result <- lapply(df, function(x) x[!is.na(x) & x != ""])
  return(result)
}

#' Detect species from gene names
#'
#' @param features Character vector of gene names
#' @return "Hs" for human, "Mm" for mouse
#' @export
detect_species <- function(features) {
  # Count different naming patterns
  ens_hs <- sum(grepl("^ENSG", features))
  ens_mm <- sum(grepl("^ENSMUSG", features))
  uppercase <- sum(features == toupper(features) & !grepl("^ENS", features))
  titlecase <- sum(features == tools::toTitleCase(tolower(features)) & 
                     !grepl("^Ens", features))
  
  # Determine species
  if (ens_hs > ens_mm && ens_hs > uppercase) {
    return("Hs")
  } else if (ens_mm > ens_hs && ens_mm > titlecase) {
    return("Mm")
  } else if (uppercase > titlecase) {
    return("Hs")
  } else {
    return("Mm")
  }
}

#' Clip values to quantile range
#'
#' @param x Numeric vector
#' @param lower_quantile Lower quantile for clipping (0-1)
#' @param upper_quantile Upper quantile for clipping (0-1)
#' @return Clipped vector
#' @export
clip_to_quantile <- function(x, lower_quantile = 0, upper_quantile = 1) {
  lower_val <- quantile(x, lower_quantile, na.rm = TRUE)
  upper_val <- quantile(x, upper_quantile, na.rm = TRUE)
  
  x[x < lower_val] <- lower_val
  x[x > upper_val] <- upper_val
  
  return(x)
}

#' Calculate Jaccard similarity between two vectors
#'
#' @param a First vector
#' @param b Second vector
#' @return Jaccard similarity (0-1)
#' @export
jaccard_similarity <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(union(a, b))
  
  if (union == 0) return(0)
  return(intersection / union)
}

#' Calculate Jaccard similarity matrix for list of gene sets
#'
#' @param gene_sets Named list of gene vectors
#' @return Jaccard similarity matrix
#' @export
jaccard_similarity_matrix <- function(gene_sets) {
  n <- length(gene_sets)
  mat <- matrix(0, n, n, dimnames = list(names(gene_sets), names(gene_sets)))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i, j] <- jaccard_similarity(gene_sets[[i]], gene_sets[[j]])
    }
  }
  
  return(mat)
}

#' Save figure with proper dimensions
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param dpi Resolution
#' @param config Configuration list (optional, for default paths)
#' @export
save_figure <- function(plot, filename, width = 7, height = 5, dpi = 300,
                        config = NULL) {
  if (!is.null(config)) {
    filepath <- file.path(config$paths$figures_dir, filename)
  } else {
    filepath <- filename
  }
  
  # Ensure directory exists
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  
  log_message(sprintf("Saved figure: %s", filepath))
}

#' Save table to CSV and RDS
#'
#' @param df Data frame to save
#' @param filename Base filename (without extension)
#' @param config Configuration list (optional, for default paths)
#' @export
save_table <- function(df, filename, config = NULL) {
  if (!is.null(config)) {
    base_path <- file.path(config$paths$tables_dir, filename)
  } else {
    base_path <- filename
  }
  
  # Ensure directory exists
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save as CSV
  csv_path <- paste0(base_path, ".csv")
  write.csv(df, csv_path, row.names = FALSE)
  
  # Save as RDS for R
  rds_path <- paste0(base_path, ".rds")
  saveRDS(df, rds_path)
  
  log_message(sprintf("Saved table: %s (.csv and .rds)", base_path))
}

#' Get unique count efficiently
#'
#' @param x Vector
#' @return Number of unique elements
#' @export
n_unique <- function(x) {
  length(unique(x))
}

#' Convert column to rownames
#'
#' @param df Data frame
#' @param col Column name to use as rownames
#' @return Data frame with rownames set
#' @export
col_to_rowname <- function(df, col) {
  rownames(df) <- df[[col]]
  df[[col]] <- NULL
  return(df)
}

#' Calculate PCA variance explained
#'
#' @param object Seurat object with PCA computed
#' @param reduction_name Name of PCA reduction
#' @return Data frame with PC variance explained
#' @export
get_pca_variance <- function(object, reduction_name = "pca") {
  stdev <- Seurat::Stdev(object, reduction = reduction_name)
  var_explained <- (stdev^2) / sum(stdev^2)
  cum_var <- cumsum(var_explained)
  
  data.frame(
    pc = seq_along(var_explained),
    var_explained = var_explained,
    cum_var = cum_var
  )
}

#' Get number of PCs for target variance explained
#'
#' @param object Seurat object with PCA computed
#' @param variance_threshold Target cumulative variance (0-1)
#' @param reduction_name Name of PCA reduction
#' @return Number of PCs
#' @export
get_npcs_for_variance <- function(object, variance_threshold = 0.9,
                                   reduction_name = "pca") {
  df_var <- get_pca_variance(object, reduction_name)
  n_pcs <- min(which(df_var$cum_var >= variance_threshold))
  
  log_message(sprintf("Using %d PCs (%.1f%% variance explained)",
                      n_pcs, df_var$cum_var[n_pcs] * 100))
  
  return(n_pcs)
}

#' Extract UMAP coordinates with metadata
#'
#' @param object Seurat object
#' @param reduction Name of reduction (default: "umap")
#' @param meta_cols Metadata columns to include (NULL for all)
#' @return Data frame with UMAP coordinates and metadata
#' @export
get_umap_coords <- function(object, reduction = "umap", meta_cols = NULL) {
  # Get UMAP embeddings
  umap_coords <- as.data.frame(Seurat::Embeddings(object, reduction = reduction))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  # Get metadata
  if (is.null(meta_cols)) {
    meta <- object@meta.data
  } else {
    meta <- object@meta.data[, meta_cols, drop = FALSE]
  }
  
  # Combine
  result <- cbind(umap_coords, meta)
  result$cell <- rownames(result)
  
  return(result)
}

#' Timer for code blocks
#'
#' @param expr Expression to time
#' @param msg Message to print with timing
#' @return Result of expression
#' @export
timed <- function(expr, msg = "Operation") {
  start <- Sys.time()
  result <- force(expr)
  elapsed <- difftime(Sys.time(), start, units = "secs")
  
  log_message(sprintf("%s completed in %.2f seconds", msg, elapsed))
  
  return(result)
}

#' Print session info for reproducibility
#'
#' @param config Configuration list (for seed info)
#' @export
print_session_info <- function(config = NULL) {
  cat("\n")
  cat("=" ,rep("=", 70), "\n", sep = "")
  cat("SESSION INFO\n")
  cat("=" ,rep("=", 70), "\n", sep = "")
  
  if (!is.null(config)) {
    cat(sprintf("Random seed: %d\n", config$reproducibility$seed))
  }
  
  cat(sprintf("Date: %s\n", Sys.time()))
  cat(sprintf("R version: %s\n", R.version.string))
  
  cat("\nKey packages:\n")
  pkgs <- c("Seurat", "sctransform", "dplyr", "ggplot2", "presto")
  for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
    }
  }
  
  cat("=" ,rep("=", 70), "\n\n", sep = "")
}
