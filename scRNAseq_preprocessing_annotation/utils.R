# =============================================================================
# Common utility functions for the OL-GBM analysis pipeline
# =============================================================================

# -----------------------------------------------------------------------------
# Package Loading
# -----------------------------------------------------------------------------

#' Load required packages, installing if necessary
#' @param packages Character vector of package names
#' @param bioc_packages Character vector of Bioconductor packages
#' @param github_packages Named character vector: c("package" = "user/repo")
load_packages <- function(packages = NULL, 
                          bioc_packages = NULL, 
                          github_packages = NULL) {
  
  # Helper to check and install CRAN packages
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing package: %s", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Install Bioconductor packages
  install_bioc_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing Bioconductor package: %s", pkg))
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Install GitHub packages
  install_github_if_missing <- function(pkg, repo) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing GitHub package: %s from %s", pkg, repo))
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github(repo)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Load CRAN packages
  if (!is.null(packages)) {
    invisible(lapply(packages, install_if_missing))
  }
  
  # Load Bioconductor packages
  if (!is.null(bioc_packages)) {
    invisible(lapply(bioc_packages, install_bioc_if_missing))
  }
  
  # Load GitHub packages
  if (!is.null(github_packages)) {
    invisible(mapply(install_github_if_missing, 
                     names(github_packages), 
                     github_packages))
  }
}

# -----------------------------------------------------------------------------
# Configuration Management
# -----------------------------------------------------------------------------

#' Load configuration from YAML file
#' @param config_path Path to config.yml
#' @return List with configuration parameters
load_config <- function(config_path = "config.yml") {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found: %s", config_path))
  }
  
  if (!requireNamespace("yaml", quietly = TRUE)) {
    install.packages("yaml")
  }
  
  config <- yaml::read_yaml(config_path)
  message("Configuration loaded successfully")
  return(config)
}

#' Create output directories from config
#' @param config Configuration list
create_directories <- function(config) {
  dirs <- c(
    config$paths$data_raw,
    config$paths$data_processed,
    config$paths$results_figures,
    config$paths$results_tables,
    config$paths$cache_dir
  )
  
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("Created directory: %s", d))
    }
  }
}

# -----------------------------------------------------------------------------
# Logging Functions
# -----------------------------------------------------------------------------

#' Print formatted message with timestamp
#' @param ... Message components
#' @param verbose Whether to print (default TRUE)
miko_message <- function(..., verbose = TRUE) {
  if (verbose) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
    message(msg)
  }
}

#' Initialize analysis log
#' @param analysis_name Name of the analysis
#' @return Data frame for logging
initiate_log <- function(analysis_name) {
  data.frame(
    timestamp = character(),
    parameter = character(),
    value = character(),
    stringsAsFactors = FALSE
  )
}

#' Add entry to analysis log
#' @param parameter Parameter name
#' @param value Parameter value
#' @param log Existing log data frame
#' @return Updated log
add_log_entry <- function(parameter, value, log) {
  new_entry <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    parameter = parameter,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  rbind(log, new_entry)
}

# -----------------------------------------------------------------------------
# Data Validation
# -----------------------------------------------------------------------------

#' Validate Seurat object has required components
#' @param object Seurat object
#' @param required_assays Character vector of required assay names
#' @param required_reductions Character vector of required reduction names
#' @param required_metadata Character vector of required metadata columns
validate_seurat <- function(object, 
                            required_assays = NULL,
                            required_reductions = NULL,
                            required_metadata = NULL) {
  
  errors <- c()
  
  # Check assays
  if (!is.null(required_assays)) {
    missing_assays <- setdiff(required_assays, names(object@assays))
    if (length(missing_assays) > 0) {
      errors <- c(errors, sprintf("Missing assays: %s", paste(missing_assays, collapse = ", ")))
    }
  }
  
  # Check reductions
  if (!is.null(required_reductions)) {
    missing_reductions <- setdiff(required_reductions, names(object@reductions))
    if (length(missing_reductions) > 0) {
      errors <- c(errors, sprintf("Missing reductions: %s", paste(missing_reductions, collapse = ", ")))
    }
  }
  
  # Check metadata
  if (!is.null(required_metadata)) {
    missing_metadata <- setdiff(required_metadata, colnames(object@meta.data))
    if (length(missing_metadata) > 0) {
      errors <- c(errors, sprintf("Missing metadata columns: %s", paste(missing_metadata, collapse = ", ")))
    }
  }
  
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"))
  }
  
  miko_message("Seurat object validation passed")
  return(TRUE)
}

#' Check if genes exist in object
#' @param object Seurat object
#' @param genes Character vector of gene names
#' @param warn_missing Whether to warn about missing genes
#' @return Character vector of genes found in object
check_genes <- function(object, genes, warn_missing = TRUE) {
  available_genes <- rownames(object)
  found <- intersect(genes, available_genes)
  missing <- setdiff(genes, available_genes)
  
  if (length(missing) > 0 && warn_missing) {
    warning(sprintf("Genes not found in object: %s", paste(head(missing, 10), collapse = ", ")))
    if (length(missing) > 10) {
      warning(sprintf("... and %d more", length(missing) - 10))
    }
  }
  
  return(found)
}

# -----------------------------------------------------------------------------
# Data Transformation Helpers
# -----------------------------------------------------------------------------

#' Clip values to quantile range
#' @param x Numeric vector
#' @param lower_quantile Lower quantile (default 0)
#' @param upper_quantile Upper quantile (default 1)
#' @return Clipped numeric vector
snip <- function(x, lower_quantile = 0, upper_quantile = 1) {
  lower_bound <- quantile(x, lower_quantile, na.rm = TRUE)
  upper_bound <- quantile(x, upper_quantile, na.rm = TRUE)
  pmin(pmax(x, lower_bound), upper_bound)
}

#' Scale vector to 0-1 range
#' @param x Numeric vector
#' @return Scaled numeric vector
scale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#' Convert long data frame to named list
#' @param df Data frame
#' @param group_by Column to group by
#' @param values Column with values
#' @return Named list
longdf_to_namedlist <- function(df, group_by, values) {
  split(df[[values]], df[[group_by]])
}

#' Convert column to rownames
#' @param df Data frame
#' @param col Column name to use as rownames
#' @return Data frame with rownames
col_to_rowname <- function(df, col) {
  rownames(df) <- df[[col]]
  df[[col]] <- NULL
  df
}

# -----------------------------------------------------------------------------
# Gene Name Conversion
# -----------------------------------------------------------------------------

#' Detect species from gene names
#' @param genes Character vector of gene names
#' @return Character: "Hs" for human, "Mm" for mouse
detect_species <- function(genes) {
  n_upper <- sum(genes == toupper(genes), na.rm = TRUE)
  n_title <- sum(genes == tools::toTitleCase(tolower(genes)), na.rm = TRUE)
  
  if (n_upper > n_title) {
    return("Hs")
  } else {
    return("Mm")
  }
}

#' Convert gene names between species
#' @param genes Character vector of gene names
#' @param to_species Target species ("Hs" or "Mm")
#' @return Character vector of converted gene names
convert_species <- function(genes, to_species) {
  if (to_species == "Hs") {
    return(toupper(genes))
  } else if (to_species == "Mm") {
    return(tools::toTitleCase(tolower(genes)))
  } else {
    stop("Unknown species: ", to_species)
  }
}

#' Check gene name representation (symbol vs ensembl)
#' @param genes Character vector of gene names
#' @return Character: "symbol" or "ensembl"
check_gene_rep <- function(genes) {
  ens_count <- sum(grepl("^ENS", genes))
  if (ens_count > length(genes) / 2) {
    return("ensembl")
  } else {
    return("symbol")
  }
}

# -----------------------------------------------------------------------------
# Statistical Helpers
# -----------------------------------------------------------------------------

#' Calculate proportion of variance explained by PCA
#' @param object Seurat object with PCA reduction
#' @param reduction_name Name of PCA reduction (default "pca")
#' @return Data frame with PC variance information
prop_var_pca <- function(object, reduction_name = "pca") {
  stdev <- object@reductions[[reduction_name]]@stdev
  var_explained <- stdev^2 / sum(stdev^2)
  cum_var <- cumsum(var_explained)
  
  data.frame(
    pc_id = seq_along(stdev),
    var_explained = var_explained,
    cum_var = cum_var
  )
}

#' Get optimal number of PCs based on variance threshold
#' @param object Seurat object with PCA
#' @param variance_threshold Cumulative variance threshold (default 0.9)
#' @return Integer number of PCs
get_optimal_pcs <- function(object, variance_threshold = 0.9) {

  df_var <- prop_var_pca(object)
  max(df_var$pc_id[df_var$cum_var < variance_threshold]) + 1
}

#' Calculate Jaccard similarity between two sets
#' @param set1 First set
#' @param set2 Second set
#' @return Numeric Jaccard similarity
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) return(0)
  intersection / union
}

#' Calculate Jaccard similarity matrix for list of sets
#' @param set_list Named list of sets
#' @return Matrix of Jaccard similarities
jaccard_similarity_matrix <- function(set_list) {
  n <- length(set_list)
  mat <- matrix(0, n, n)
  rownames(mat) <- colnames(mat) <- names(set_list)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i, j] <- jaccard_similarity(set_list[[i]], set_list[[j]])
    }
  }
  mat
}

# -----------------------------------------------------------------------------
# Parallel Processing Helpers
# -----------------------------------------------------------------------------

#' Set up parallel processing
#' @param n_cores Number of cores to use
#' @param max_memory Maximum memory per worker in GB
setup_parallel <- function(n_cores = 4, max_memory = 8) {
  available_cores <- parallel::detectCores()
  n_cores <- min(n_cores, available_cores)
  
  options(future.globals.maxSize = max_memory * 1024^3)
  
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = n_cores)
  }
  
  miko_message(sprintf("Parallel processing enabled with %d cores", n_cores))
  return(n_cores)
}

# -----------------------------------------------------------------------------
# Output Helpers
# -----------------------------------------------------------------------------

#' Save figure with consistent naming
#' @param plot ggplot object
#' @param name Figure name 
#' @param config Configuration list
#' @param width Figure width in inches
#' @param height Figure height in inches
save_figure <- function(plot, name, config, width = 7, height = 5) {
  formats <- config$visualization$figure_format
  dpi <- config$visualization$dpi
  base_path <- config$paths$results_figures
  
  for (fmt in formats) {
    filename <- file.path(base_path, paste0(name, ".", fmt))
    ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi)
    miko_message(sprintf("Saved figure: %s", filename))
  }
}

#' Save table with consistent naming
#' @param df Data frame to save
#' @param name Table name 
#' @param config Configuration list
#' @param formats Output formats (default c("csv", "rds"))
save_table <- function(df, name, config, formats = c("csv", "rds")) {
  base_path <- config$paths$results_tables
  
  if ("csv" %in% formats) {
    csv_path <- file.path(base_path, paste0(name, ".csv"))
    write.csv(df, csv_path, row.names = FALSE)
    miko_message(sprintf("Saved table: %s", csv_path))
  }
  
  if ("rds" %in% formats) {
    rds_path <- file.path(base_path, paste0(name, ".rds"))
    saveRDS(df, rds_path)
    miko_message(sprintf("Saved table: %s", rds_path))
  }
}

#' Save Seurat object
#' @param object Seurat object
#' @param name Object name
#' @param config Configuration list
save_seurat <- function(object, name, config) {
  rds_path <- file.path(config$paths$data_processed, paste0(name, ".rds"))
  saveRDS(object, rds_path)
  miko_message(sprintf("Saved Seurat object: %s", rds_path))
}

# -----------------------------------------------------------------------------
# Session Info
# -----------------------------------------------------------------------------

#' Print and save session info
#' @param output_path Path to save session info (optional)
print_session_info <- function(output_path = NULL) {
  si <- sessionInfo()
  print(si)
  
  if (!is.null(output_path)) {
    writeLines(capture.output(print(si)), output_path)
    miko_message(sprintf("Session info saved: %s", output_path))
  }
  
  invisible(si)
}
