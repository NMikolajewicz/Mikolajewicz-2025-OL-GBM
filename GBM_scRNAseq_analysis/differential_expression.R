# ==============================================================================
# differential_expression.R - DEG and Abundance Analysis Functions
# ==============================================================================
# Functions for differential expression analysis and cell type abundance.
# Uses presto::wilcoxauc for fast Wilcoxon tests (as per manuscript).
# ==============================================================================

#' Run differential expression using presto
#'
#' @param object Seurat object
#' @param group_by Metadata column for grouping (e.g., "seurat_clusters")
#' @param config Configuration list
#' @param assay Assay to use
#' @return Data frame with DE results
#' @export
run_de_presto <- function(object, group_by = "seurat_clusters", 
                           config = NULL, assay = NULL) {
  
  log_message(sprintf("Running differential expression (presto) by %s...", 
                      group_by))
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Get expression matrix
  expr_mat <- Seurat::GetAssayData(object, slot = "data", assay = assay)
  
  # Get group labels
  groups <- object@meta.data[[group_by]]
  
  # Run wilcoxauc
  de_results <- presto::wilcoxauc(
    expr_mat,
    groups
  )
  
  # Add adjusted p-values
  de_results$padj <- p.adjust(de_results$pval, method = "BH")
  
  # Filter by thresholds if config provided
  if (!is.null(config)) {
    de_results$significant <- 
      de_results$padj < config$differential_expression$fdr_threshold &
      abs(de_results$logFC) >= config$differential_expression$logfc_threshold &
      de_results$pct_in >= config$differential_expression$min_pct * 100
  }
  
  log_message(sprintf("  Found %d significant DE genes", 
                      sum(de_results$significant, na.rm = TRUE)))
  
  return(de_results)
}

#' Get top markers per cluster
#'
#' @param de_results DE results from run_de_presto
#' @param n_top Number of top markers per group
#' @param direction "up", "down", or "both"
#' @return Data frame with top markers
#' @export
get_top_markers <- function(de_results, n_top = 20, direction = "up") {
  
  if (direction == "up") {
    de_filtered <- de_results[de_results$logFC > 0, ]
  } else if (direction == "down") {
    de_filtered <- de_results[de_results$logFC < 0, ]
  } else {
    de_filtered <- de_results
  }
  
  # Get top markers per group
  top_markers <- de_filtered %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(padj, dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup()
  
  return(top_markers)
}

#' Run pairwise differential expression
#'
#' @param object Seurat object
#' @param group1 Cells in group 1 (column value or cell names)
#' @param group2 Cells in group 2 (column value or cell names)
#' @param group_col Column name if using column values
#' @param assay Assay to use
#' @return Data frame with DE results
#' @export
run_de_pairwise <- function(object, group1, group2, 
                             group_col = NULL, assay = NULL) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Get cells for each group
  if (!is.null(group_col)) {
    cells1 <- colnames(object)[object@meta.data[[group_col]] == group1]
    cells2 <- colnames(object)[object@meta.data[[group_col]] == group2]
  } else {
    cells1 <- group1
    cells2 <- group2
  }
  
  # Get expression matrix
  expr_mat <- Seurat::GetAssayData(object, slot = "data", assay = assay)
  
  # Create group labels
  all_cells <- c(cells1, cells2)
  labels <- c(rep("group1", length(cells1)), rep("group2", length(cells2)))
  
  # Subset and run
  de_results <- presto::wilcoxauc(
    expr_mat[, all_cells],
    labels
  )
  
  # Filter to group1 (test group)
  de_results <- de_results[de_results$group == "group1", ]
  de_results$padj <- p.adjust(de_results$pval, method = "BH")
  
  return(de_results)
}

#' Calculate cell type fractional abundance
#'
#' @param object Seurat object
#' @param cell_type_col Column with cell type annotations
#' @param sample_col Column with sample identifiers
#' @return Data frame with fractional abundances
#' @export
calculate_abundance <- function(object, cell_type_col = "cell_type",
                                 sample_col = "sample") {
  
  log_message("Calculating cell type abundances...")
  
  df_meta <- object@meta.data[, c(cell_type_col, sample_col)]
  colnames(df_meta) <- c("cell_type", "sample")
  
  # Count cells per type per sample
  df_counts <- df_meta %>%
    dplyr::group_by(sample, cell_type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  # Calculate totals per sample
  df_totals <- df_meta %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(total = dplyr::n(), .groups = "drop")
  
  # Calculate fractions
  df_abundance <- dplyr::left_join(df_counts, df_totals, by = "sample")
  df_abundance$fraction <- df_abundance$n / df_abundance$total
  
  return(df_abundance)
}

#' Compare abundances between conditions (paired or unpaired)
#'
#' @param df_abundance Abundance data frame from calculate_abundance
#' @param condition_col Column with condition info
#' @param conditions Vector of two conditions to compare
#' @param pair_col Column for pairing (NULL for unpaired)
#' @return Data frame with comparison results
#' @export
compare_abundance <- function(df_abundance, condition_col, conditions,
                               pair_col = NULL) {
  
  log_message("Comparing abundances between conditions...")
  
  # Get condition info
  df_abundance$condition <- df_abundance[[condition_col]]
  
  # Filter to conditions of interest
  df_filt <- df_abundance[df_abundance$condition %in% conditions, ]
  
  cell_types <- unique(df_filt$cell_type)
  results <- list()
  
  for (ct in cell_types) {
    df_ct <- df_filt[df_filt$cell_type == ct, ]
    
    group1 <- df_ct$fraction[df_ct$condition == conditions[1]]
    group2 <- df_ct$fraction[df_ct$condition == conditions[2]]
    
    if (!is.null(pair_col)) {
      # Paired test
      df_wide <- df_ct %>%
        tidyr::pivot_wider(
          id_cols = all_of(pair_col),
          names_from = condition,
          values_from = fraction
        )
      
      if (all(conditions %in% colnames(df_wide))) {
        test_result <- wilcox.test(
          df_wide[[conditions[1]]],
          df_wide[[conditions[2]]],
          paired = TRUE
        )
        p_paired <- test_result$p.value
      } else {
        p_paired <- NA
      }
    } else {
      p_paired <- NA
    }
    
    # Unpaired test
    test_unpaired <- wilcox.test(group1, group2, paired = FALSE)
    
    results[[ct]] <- data.frame(
      cell_type = ct,
      mean_1 = mean(group1, na.rm = TRUE),
      mean_2 = mean(group2, na.rm = TRUE),
      log2fc = log2(mean(group2, na.rm = TRUE) / mean(group1, na.rm = TRUE)),
      p_unpaired = test_unpaired$p.value,
      p_paired = p_paired
    )
  }
  
  df_results <- dplyr::bind_rows(results)
  df_results$padj_unpaired <- p.adjust(df_results$p_unpaired, method = "BH")
  df_results$padj_paired <- p.adjust(df_results$p_paired, method = "BH")
  
  return(df_results)
}

#' Plot fractional abundance comparison
#'
#' @param df_abundance Abundance data frame
#' @param cell_type Cell type to plot
#' @param condition_col Column with condition info
#' @param pair_col Column for paired lines (optional)
#' @return ggplot object
#' @export
plot_abundance_comparison <- function(df_abundance, cell_type,
                                       condition_col, pair_col = NULL) {
  
  df_plot <- df_abundance[df_abundance$cell_type == cell_type, ]
  df_plot$condition <- df_plot[[condition_col]]
  
  p <- ggplot2::ggplot(df_plot, 
                       ggplot2::aes(x = condition, y = fraction)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition), alpha = 0.5) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_y_log10() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = condition_col,
      y = "Fractional Abundance",
      title = cell_type
    ) +
    ggplot2::theme(legend.position = "none")
  
  if (!is.null(pair_col)) {
    df_plot$pair <- df_plot[[pair_col]]
    p <- p + ggplot2::geom_line(ggplot2::aes(group = pair), alpha = 0.5)
  }
  
  return(p)
}

#' Run module score enrichment
#'
#' @param object Seurat object
#' @param gene_lists Named list of gene vectors
#' @param scale Whether to scale scores
#' @param n_ctrl Number of control genes
#' @return Data frame with module scores
#' @export
run_module_scoring <- function(object, gene_lists, scale = TRUE, n_ctrl = 100) {
  
  log_message(sprintf("Calculating module scores for %d gene sets...", 
                      length(gene_lists)))
  
  # Filter to available genes
  gene_lists_filt <- lapply(gene_lists, function(x) {
    intersect(x, rownames(object))
  })
  
  # Remove empty sets
  gene_lists_filt <- gene_lists_filt[sapply(gene_lists_filt, length) >= 5]
  
  # Add module scores
  for (name in names(gene_lists_filt)) {
    object <- Seurat::AddModuleScore(
      object,
      features = list(gene_lists_filt[[name]]),
      name = name,
      ctrl = n_ctrl
    )
  }
  
  # Extract scores
  score_cols <- paste0(names(gene_lists_filt), "1")
  df_scores <- object@meta.data[, score_cols, drop = FALSE]
  colnames(df_scores) <- names(gene_lists_filt)
  
  if (scale) {
    df_scores <- as.data.frame(scale(df_scores))
  }
  
  return(df_scores)
}

#' Aggregate expression by cluster
#'
#' @param object Seurat object
#' @param group_by Column for grouping
#' @param assay Assay to use
#' @param slot Slot to use (data, counts, scale.data)
#' @param method Aggregation method (mean, median, sum)
#' @return Aggregated expression matrix (genes x clusters)
#' @export
aggregate_expression <- function(object, group_by = "seurat_clusters",
                                   assay = NULL, slot = "data",
                                   method = "mean") {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  expr_mat <- Seurat::GetAssayData(object, slot = slot, assay = assay)
  groups <- object@meta.data[[group_by]]
  unique_groups <- sort(unique(groups))
  
  agg_mat <- matrix(
    NA,
    nrow = nrow(expr_mat),
    ncol = length(unique_groups),
    dimnames = list(rownames(expr_mat), unique_groups)
  )
  
  for (g in unique_groups) {
    cells <- which(groups == g)
    if (method == "mean") {
      agg_mat[, g] <- Matrix::rowMeans(expr_mat[, cells, drop = FALSE])
    } else if (method == "median") {
      agg_mat[, g] <- apply(expr_mat[, cells, drop = FALSE], 1, median)
    } else if (method == "sum") {
      agg_mat[, g] <- Matrix::rowSums(expr_mat[, cells, drop = FALSE])
    }
  }
  
  return(agg_mat)
}

#' Create dot plot data for markers
#'
#' @param object Seurat object
#' @param features Genes to plot
#' @param group_by Grouping variable
#' @param assay Assay to use
#' @return Data frame for dot plot
#' @export
prepare_dotplot_data <- function(object, features, group_by = "seurat_clusters",
                                  assay = NULL) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  expr_mat <- Seurat::GetAssayData(object, slot = "data", assay = assay)
  
  # Filter to available features
  features <- intersect(features, rownames(expr_mat))
  expr_mat <- expr_mat[features, , drop = FALSE]
  
  groups <- object@meta.data[[group_by]]
  unique_groups <- sort(unique(groups))
  
  # Calculate mean expression and percent expressing
  df_list <- list()
  
  for (g in unique_groups) {
    cells <- which(groups == g)
    
    avg_exp <- Matrix::rowMeans(expr_mat[, cells, drop = FALSE])
    pct_exp <- Matrix::rowMeans(expr_mat[, cells, drop = FALSE] > 0) * 100
    
    df_list[[as.character(g)]] <- data.frame(
      feature = features,
      group = g,
      avg.exp = avg_exp,
      pct.exp = pct_exp
    )
  }
  
  df_dot <- dplyr::bind_rows(df_list)
  
  # Scale average expression within features
  df_dot <- df_dot %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(avg.exp.scaled = scale(avg.exp)[, 1]) %>%
    dplyr::ungroup()
  
  return(df_dot)
}

#' Create marker dot plot
#'
#' @param df_dot Dot plot data from prepare_dotplot_data
#' @param feature_order Order of features (top to bottom)
#' @param group_order Order of groups (left to right)
#' @return ggplot object
#' @export
plot_markers_dot <- function(df_dot, feature_order = NULL, group_order = NULL) {
  
  if (!is.null(feature_order)) {
    df_dot$feature <- factor(df_dot$feature, levels = rev(feature_order))
  }
  
  if (!is.null(group_order)) {
    df_dot$group <- factor(df_dot$group, levels = group_order)
  }
  
  # Get color limits
  color_lim <- max(abs(df_dot$avg.exp.scaled), na.rm = TRUE)
  
  p <- ggplot2::ggplot(df_dot, 
                       ggplot2::aes(x = group, y = feature)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::scale_color_distiller(
      palette = "RdBu",
      limits = c(-color_lim, color_lim)
    ) +
    ggplot2::scale_size_continuous(range = c(0, 5)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Cluster",
      y = "Gene",
      color = "Scaled\nExpression",
      size = "% Expressing"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_line(color = "grey90")
    )
  
  return(p)
}
