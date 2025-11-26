# ==============================================================================
# cnv.R - Copy Number Variation Analysis for Tumor Calling
# ==============================================================================

#' Get chromosome position gene sets from MSigDB C1
#'
#' @param species "Homo sapiens" or "Mus musculus"
#' @param exclude_hla Exclude HLA genes from chr6 (recommended)
#' @return Named list of gene sets by chromosome arm
#' @export
get_chromosome_gene_sets <- function(species = "Homo sapiens",
                                      exclude_hla = TRUE) {
  
  log_message("Loading chromosome position gene sets...")
  
  # Get C1 (positional) gene sets
  chr_pos <- msigdbr::msigdbr(species = species, category = "C1")
  
  # Filter HLA genes if requested (per manuscript: HLA can manifest as CNV)
  if (exclude_hla) {
    chr_pos <- chr_pos[!grepl("^HLA-", chr_pos$gene_symbol), ]
    log_message("  Excluded HLA genes from chromosome 6")
  }
  
  # Convert to named list
  chr_list <- longDF_to_named_list(chr_pos, 
                                    group_by = "gs_name", 
                                    values = "gene_symbol")
  
  log_message(sprintf("  Loaded %d chromosome bands", length(chr_list)))
  
  return(chr_list)
}

#' Aggregate chromosome band gene sets to chromosome arms
#'
#' @param chr_band_list Named list of chromosome band gene sets
#' @return Named list of chromosome arm gene sets (e.g., "1p", "1q", ...)
#' @export
aggregate_to_chromosome_arms <- function(chr_band_list) {
  
  # Parse chromosome info from band names
  df_chr <- data.frame(
    band = names(chr_band_list)
  )
  
  df_chr$chr_num <- gsub("chr", "", stringr::str_extract(df_chr$band, "chr[0-9XYM]*"))
  df_chr$chr_arm <- stringr::str_extract(df_chr$band, "[pq]")
  df_chr$chr_num_arm <- paste0(df_chr$chr_num, df_chr$chr_arm)
  
  # Filter valid entries
  df_chr <- df_chr[complete.cases(df_chr), ]
  df_chr <- df_chr[df_chr$chr_num != "Y", ]  # Exclude Y chromosome
  
  # Aggregate genes by chromosome arm
  chr_arms <- unique(df_chr$chr_num_arm)
  chr_arm_list <- list()
  
  for (arm in chr_arms) {
    bands <- df_chr$band[df_chr$chr_num_arm == arm]
    genes <- unique(unlist(chr_band_list[bands]))
    chr_arm_list[[arm]] <- genes
  }
  
  # Order chromosomes properly
  chr_order <- c(paste0(rep(1:22, each = 2), rep(c("p", "q"), 22)), "Xp", "Xq")
  chr_arm_list <- chr_arm_list[intersect(chr_order, names(chr_arm_list))]
  
  log_message(sprintf("  Aggregated to %d chromosome arms", length(chr_arm_list)))
  
  return(chr_arm_list)
}

#' Calculate chromosome arm expression scores
#'
#' @param object Seurat object
#' @param chr_arm_list Named list of chromosome arm gene sets
#' @param scale Whether to scale scores
#' @param n_ctrl Number of control genes for module scoring
#' @return Data frame with cells x chromosome arms
#' @export
score_chromosome_arms <- function(object, chr_arm_list, 
                                   scale = TRUE, n_ctrl = 100) {
  
  log_message("Scoring chromosome arm expression...")
  
  # Filter gene sets to available genes
  available_genes <- rownames(object)
  chr_arm_list_filt <- lapply(chr_arm_list, function(x) {
    intersect(x, available_genes)
  })
  
  # Remove empty gene sets
  chr_arm_list_filt <- chr_arm_list_filt[
    sapply(chr_arm_list_filt, length) >= 10
  ]
  
  log_message(sprintf("  Using %d chromosome arms with sufficient genes",
                      length(chr_arm_list_filt)))
  
  # Score each chromosome arm
  object_copy <- object
  for (arm_name in names(chr_arm_list_filt)) {
    object_copy <- Seurat::AddModuleScore(
      object_copy,
      features = list(chr_arm_list_filt[[arm_name]]),
      name = arm_name,
      ctrl = n_ctrl
    )
  }
  
  # Extract scores
  score_cols <- paste0(names(chr_arm_list_filt), "1")
  df_scores <- object_copy@meta.data[, score_cols, drop = FALSE]
  colnames(df_scores) <- names(chr_arm_list_filt)
  
  # Scale if requested
  if (scale) {
    df_scores <- as.data.frame(scale(df_scores))
  }
  
  return(df_scores)
}

#' Calculate pairwise chromosome arm ratios
#'
#' @param df_scores Data frame of chromosome arm scores (cells x arms)
#' @return Data frame of pairwise ratios
#' @export
calculate_chr_ratios <- function(df_scores) {
  
  log_message("Calculating pairwise chromosome ratios...")
  
  arms <- colnames(df_scores)
  n_arms <- length(arms)
  
  # Initialize result matrix
  df_ratios <- data.frame(matrix(NA, nrow = nrow(df_scores), ncol = 0))
  rownames(df_ratios) <- rownames(df_scores)
  
  # Calculate log ratios between each pair
  for (i in seq_len(n_arms)) {
    for (j in seq_len(n_arms)) {
      if (i == j) next
      
      ratio_name <- paste0("c", arms[i], "_c", arms[j])
      
      # Log ratio (add small constant to avoid log(0))
      ratio <- log((exp(df_scores[, arms[i]]) + 0.01) / 
                   (exp(df_scores[, arms[j]]) + 0.01))
      
      df_ratios[[ratio_name]] <- ratio
    }
  }
  
  log_message(sprintf("  Computed %d pairwise ratios", ncol(df_ratios)))
  
  return(df_ratios)
}

#' Correct CNV scores using normal reference cells
#'
#' @param df_scores Data frame of CNV scores (cells x features)
#' @param normal_cells Logical or character vector of normal reference cells
#' @return Corrected data frame
#' @export
correct_cnv_scores <- function(df_scores, normal_cells) {
  
  log_message("Correcting CNV scores using normal reference...")
  
  # Get normal cell indices
  if (is.logical(normal_cells)) {
    normal_idx <- which(normal_cells)
  } else {
    normal_idx <- which(rownames(df_scores) %in% normal_cells)
  }
  
  log_message(sprintf("  Using %d normal reference cells", length(normal_idx)))
  
  # Calculate mean scores in normal cells
  normal_means <- colMeans(df_scores[normal_idx, , drop = FALSE], na.rm = TRUE)
  
  # Subtract normal mean from all cells
  df_corrected <- sweep(df_scores, 2, normal_means, "-")
  
  return(df_corrected)
}

#' Identify microglia cells for CNV reference
#'
#' @param object Seurat object
#' @param microglia_clusters Cluster IDs identified as microglia
#' @param microglia_markers Marker genes for microglia
#' @return Logical vector indicating microglia cells
#' @export
identify_microglia <- function(object, microglia_clusters = NULL,
                                microglia_markers = c("P2RY12", "CD163")) {
  
  if (!is.null(microglia_clusters)) {
    # Use cluster-based identification
    is_microglia <- Seurat::Idents(object) %in% microglia_clusters
    
  } else {
    # Use marker-based identification
    markers_present <- intersect(microglia_markers, rownames(object))
    
    if (length(markers_present) == 0) {
      stop("No microglia markers found in object")
    }
    
    # Score cells for microglia markers
    object <- Seurat::AddModuleScore(
      object,
      features = list(markers_present),
      name = "microglia_score"
    )
    
    # Threshold at mean + SD
    scores <- object$microglia_score1
    threshold <- mean(scores) + sd(scores)
    is_microglia <- scores > threshold
  }
  
  log_message(sprintf("Identified %d microglia cells (%.1f%%)",
                      sum(is_microglia),
                      sum(is_microglia) / ncol(object) * 100))
  
  return(is_microglia)
}

#' Run PCA on CNV matrix for tumor calling
#'
#' @param df_cnv Corrected CNV score matrix
#' @param n_pcs Number of PCs to compute
#' @return List with embeddings, loadings, and variance explained
#' @export
run_cnv_pca <- function(df_cnv, n_pcs = 10) {
  
  log_message("Running PCA on CNV matrix...")
  
  # Remove NA/Inf values
  df_cnv <- df_cnv[complete.cases(df_cnv), ]
  df_cnv <- df_cnv[apply(df_cnv, 1, function(x) all(is.finite(x))), ]
  
  # Run PCA
  pca_result <- prcomp(df_cnv, scale. = TRUE, center = TRUE)
  
  # Get variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  
  # Get top loadings for PC1 and PC2
  loadings_pc1 <- sort(abs(pca_result$rotation[, 1]), decreasing = TRUE)[1:5]
  loadings_pc2 <- sort(abs(pca_result$rotation[, 2]), decreasing = TRUE)[1:5]
  
  log_message(sprintf("  PC1: %.1f%% variance", var_explained[1] * 100))
  log_message(sprintf("  PC2: %.1f%% variance", var_explained[2] * 100))
  log_message(sprintf("  Top PC1 features: %s", 
                      paste(names(loadings_pc1), collapse = ", ")))
  
  return(list(
    embeddings = pca_result$x[, 1:min(n_pcs, ncol(pca_result$x))],
    loadings = pca_result$rotation[, 1:min(n_pcs, ncol(pca_result$rotation))],
    var_explained = var_explained[1:min(n_pcs, length(var_explained))],
    cells = rownames(df_cnv)
  ))
}

#' Score cells for tumor markers
#'
#' @param object Seurat object
#' @param tumor_markers Vector of tumor marker gene names
#' @param scale Whether to scale scores
#' @return Numeric vector of tumor scores
#' @export
score_tumor_markers <- function(object, tumor_markers, scale = TRUE) {
  
  log_message("Scoring tumor markers...")
  
  # Filter to available markers
  markers_present <- intersect(tumor_markers, rownames(object))
  
  log_message(sprintf("  Using %d/%d markers", 
                      length(markers_present), length(tumor_markers)))
  
  if (length(markers_present) < 3) {
    warning("Less than 3 tumor markers found")
    return(rep(0, ncol(object)))
  }
  
  # Calculate eigengene (PC1 of marker expression)
  expr_mat <- Seurat::GetAssayData(object, slot = "data")[markers_present, ]
  expr_mat <- as.matrix(t(expr_mat))
  
  # Scale
  expr_mat <- scale(expr_mat)
  
  # Get first PC
  pca <- prcomp(expr_mat, scale. = FALSE, center = FALSE)
  scores <- pca$x[, 1]
  
  # Align direction (higher = more tumor-like)
  mean_expr <- rowMeans(expr_mat)
  if (cor(scores, mean_expr) < 0) {
    scores <- -scores
  }
  
  if (scale) {
    scores <- scale(scores)[, 1]
  }
  
  return(scores)
}

#' Classify cells as tumor or normal using GMM
#'
#' @param scores Numeric vector of scores (CNV or marker-based)
#' @param n_clusters Number of clusters for GMM (default 2)
#' @return List with classifications and probabilities
#' @export
classify_cells_gmm <- function(scores, n_clusters = 2) {
  
  # Fit GMM
  gmm <- mclust::Mclust(scores, G = seq_len(n_clusters), verbose = FALSE)
  
  # Get classifications
  classifications <- gmm$classification
  
  # Get uncertainty
  uncertainty <- gmm$uncertainty
  
  # Get cluster means to determine which is "tumor"
  cluster_means <- tapply(scores, classifications, mean)
  tumor_cluster <- which.max(cluster_means)
  
  # Calculate probability of being tumor
  if (!is.null(gmm$z)) {
    prob_tumor <- gmm$z[, tumor_cluster]
  } else {
    prob_tumor <- as.numeric(classifications == tumor_cluster)
  }
  
  return(list(
    classification = classifications,
    uncertainty = uncertainty,
    prob_tumor = prob_tumor,
    tumor_cluster = tumor_cluster
  ))
}

#' Calculate posterior probability of being tumor
#'
#' @param cnv_prob Probability from CNV analysis
#' @param marker_prob Probability from marker analysis
#' @return Posterior probability (1 - product of complements)
#' @export
calculate_posterior_prob <- function(cnv_prob, marker_prob) {
  
  # Bayesian combination: P(tumor) = 1 - (1-p1)(1-p2)
  posterior <- 1 - (1 - cnv_prob) * (1 - marker_prob)
  
  return(posterior)
}

#' Complete CNV-based tumor calling pipeline
#'
#' @param object Seurat object
#' @param config Configuration list
#' @param microglia_clusters Cluster IDs for microglia (NULL for auto-detect)
#' @return Seurat object with tumor annotations
#' @export
call_tumor_cells <- function(object, config, microglia_clusters = NULL) {
  
  log_message("\n=== CNV-based Tumor Calling ===")
  
  # Use reference clusters from config if not provided
  if (is.null(microglia_clusters)) {
    microglia_clusters <- config$cnv$reference_clusters
  }
  
  # 1. Identify microglia (normal reference)
  is_microglia <- identify_microglia(object, microglia_clusters)
  object$is_microglia <- is_microglia
  
  # 2. Get chromosome gene sets
  chr_band_list <- get_chromosome_gene_sets()
  chr_arm_list <- aggregate_to_chromosome_arms(chr_band_list)
  
  # 3. Score chromosome arms
  df_chr_scores <- score_chromosome_arms(object, chr_arm_list, scale = TRUE)
  
  # 4. Correct using microglia
  df_chr_corrected <- correct_cnv_scores(df_chr_scores, is_microglia)
  
  # 5. Run PCA on CNV matrix
  cnv_pca <- run_cnv_pca(df_chr_corrected)
  
  # 6. Classify using GMM on PC1+PC2
  cnv_scores <- sqrt(cnv_pca$embeddings[, 1]^2 + cnv_pca$embeddings[, 2]^2)
  names(cnv_scores) <- cnv_pca$cells
  
  cnv_class <- classify_cells_gmm(cnv_scores)
  
  # 7. Score tumor markers
  tumor_markers <- config$tumor_markers$gsc_markers
  marker_scores <- score_tumor_markers(object, tumor_markers)
  marker_class <- classify_cells_gmm(marker_scores)
  
  # 8. Calculate posterior probability
  # Align probabilities to common cells
  common_cells <- intersect(names(cnv_class$prob_tumor), colnames(object))
  
  cnv_prob <- rep(0, ncol(object))
  names(cnv_prob) <- colnames(object)
  cnv_prob[common_cells] <- cnv_class$prob_tumor[common_cells]
  
  posterior_prob <- calculate_posterior_prob(cnv_prob, marker_class$prob_tumor)
  
  # 9. Add to metadata
  threshold <- config$cnv$posterior_prob_threshold
  
  object$cnv_prob <- cnv_prob
  object$marker_prob <- marker_class$prob_tumor
  object$tumor_prob <- posterior_prob
  object$is_tumor <- posterior_prob > threshold
  
  # Summary
  n_tumor <- sum(object$is_tumor)
  n_normal <- sum(!object$is_tumor)
  
  log_message(sprintf("\nTumor calling results:"))
  log_message(sprintf("  Tumor cells: %d (%.1f%%)", 
                      n_tumor, n_tumor / ncol(object) * 100))
  log_message(sprintf("  Normal cells: %d (%.1f%%)",
                      n_normal, n_normal / ncol(object) * 100))
  
  return(object)
}

#' Plot CNV heatmap
#'
#' @param df_chr_scores Chromosome arm scores matrix
#' @param cell_order Optional ordering of cells
#' @param chr_order Ordering of chromosomes
#' @return ggplot object
#' @export
plot_cnv_heatmap <- function(df_chr_scores, cell_order = NULL, 
                              chr_order = NULL) {
  
  # Prepare data
  df_long <- tidyr::pivot_longer(
    tibble::rownames_to_column(df_chr_scores, "cell"),
    cols = -cell,
    names_to = "chr_arm",
    values_to = "score"
  )
  
  # Order cells
  if (!is.null(cell_order)) {
    df_long$cell <- factor(df_long$cell, levels = cell_order)
  }
  
  # Order chromosomes
  if (!is.null(chr_order)) {
    df_long$chr_arm <- factor(df_long$chr_arm, levels = chr_order)
  }
  
  # Clip values for visualization
  df_long$score_clipped <- clip_to_quantile(df_long$score, 0.01, 0.99)
  
  # Plot
  ggplot2::ggplot(df_long, ggplot2::aes(x = chr_arm, y = cell, fill = score_clipped)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      name = "CNV Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)
    ) +
    ggplot2::labs(x = "Chromosome Arm", y = "Cells")
}
