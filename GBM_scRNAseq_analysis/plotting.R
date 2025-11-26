# ==============================================================================
# plotting.R - Visualization Functions 
# ==============================================================================

#' Custom theme for publication figures
#'
#' @param base_size Base font size
#' @param legend Show legend
#' @param grid Show grid
#' @param x_axis_rotation X-axis label rotation
#' @return ggplot2 theme
#' @export
theme_publication <- function(base_size = 11, legend = TRUE, 
                               grid = FALSE, x_axis_rotation = 0) {
  
  theme_base <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # Text
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      
      # Axis
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      
      # Legend
      legend.position = if (legend) "right" else "none",
      legend.background = ggplot2::element_blank(),
      
      # Panel
      panel.border = ggplot2::element_blank()
    )
  
  if (!grid) {
    theme_base <- theme_base + 
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  if (x_axis_rotation != 0) {
    theme_base <- theme_base +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = x_axis_rotation, 
          hjust = 1, 
          vjust = 1
        )
      )
  }
  
  return(theme_base)
}

#' Create UMAP plot with cell annotations
#'
#' @param object Seurat object
#' @param group_by Metadata column for coloring
#' @param reduction Reduction to use
#' @param pt_size Point size
#' @param label Show cluster labels
#' @param colors Custom color palette (optional)
#' @param title Plot title
#' @return ggplot object
#' @export
plot_umap <- function(object, group_by = "seurat_clusters", 
                       reduction = "umap", pt_size = 0.5,
                       label = TRUE, colors = NULL, title = NULL) {
  
  # Get coordinates
  coords <- get_umap_coords(object, reduction = reduction)
  coords$group <- object@meta.data[[group_by]]
  
  # Create plot
  p <- ggplot2::ggplot(coords, 
                       ggplot2::aes(x = UMAP_1, y = UMAP_2, color = group)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.7) +
    theme_publication() +
    ggplot2::labs(
      x = "UMAP 1",
      y = "UMAP 2",
      color = group_by,
      title = title
    ) +
    ggplot2::coord_fixed()
  
  # Add custom colors
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }
  
  # Add labels
  if (label) {
    label_coords <- coords %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2),
        .groups = "drop"
      )
    
    p <- p + ggplot2::geom_text(
      data = label_coords,
      ggplot2::aes(label = group),
      color = "black",
      size = 3,
      fontface = "bold"
    )
  }
  
  return(p)
}

#' Create feature expression UMAP
#'
#' @param object Seurat object
#' @param feature Gene name or metadata column
#' @param reduction Reduction to use
#' @param pt_size Point size
#' @param order Order points by expression
#' @return ggplot object
#' @export
plot_feature_umap <- function(object, feature, reduction = "umap",
                               pt_size = 0.5, order = TRUE) {
  
  # Get coordinates
  coords <- get_umap_coords(object, reduction = reduction)
  
  # Get feature values
  if (feature %in% rownames(object)) {
    coords$value <- Seurat::FetchData(object, vars = feature)[, 1]
  } else if (feature %in% colnames(object@meta.data)) {
    coords$value <- object@meta.data[[feature]]
  } else {
    stop(sprintf("Feature '%s' not found", feature))
  }
  
  # Order by expression
  if (order) {
    coords <- coords[order(coords$value), ]
  }
  
  p <- ggplot2::ggplot(coords, 
                       ggplot2::aes(x = UMAP_1, y = UMAP_2, color = value)) +
    ggplot2::geom_point(size = pt_size) +
    viridis::scale_color_viridis(name = feature) +
    theme_publication() +
    ggplot2::labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = feature
    ) +
    ggplot2::coord_fixed()
  
  return(p)
}

#' Create heatmap with pheatmap
#'
#' @param mat Matrix to plot (features x samples)
#' @param scale Scale rows, columns, or none
#' @param cluster_rows Cluster rows
#' @param cluster_cols Cluster columns
#' @param annotation_col Column annotations (data.frame)
#' @param color Color palette
#' @param ... Additional arguments to pheatmap
#' @return pheatmap object
#' @export
plot_heatmap <- function(mat, scale = "row", 
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          annotation_col = NULL, color = NULL, ...) {
  
  if (is.null(color)) {
    color <- colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    )(100)
  }
  
  pheatmap::pheatmap(
    mat,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_col = annotation_col,
    color = color,
    border_color = NA,
    ...
  )
}

#' Create volcano plot for DE results
#'
#' @param de_results DE results data frame
#' @param logfc_col Column with log fold change
#' @param pval_col Column with p-values or adjusted p-values
#' @param logfc_thresh Fold change threshold for significance
#' @param pval_thresh P-value threshold for significance
#' @param label_genes Genes to label (optional)
#' @param n_top_label Number of top genes to label
#' @return ggplot object
#' @export
plot_volcano <- function(de_results, 
                          logfc_col = "logFC", 
                          pval_col = "padj",
                          logfc_thresh = 0.5,
                          pval_thresh = 0.05,
                          label_genes = NULL,
                          n_top_label = 10) {
  
  # Prepare data
  df <- data.frame(
    gene = de_results$feature,
    logfc = de_results[[logfc_col]],
    pval = -log10(de_results[[pval_col]] + 1e-300)
  )
  
  # Define significance
  df$significance <- "NS"
  df$significance[df$pval > -log10(pval_thresh) & df$logfc > logfc_thresh] <- "Up"
  df$significance[df$pval > -log10(pval_thresh) & df$logfc < -logfc_thresh] <- "Down"
  
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = logfc, y = pval, color = significance)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::scale_color_manual(
      values = c("Up" = "red", "Down" = "blue", "NS" = "grey")
    ) +
    ggplot2::geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), 
                        linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(pval_thresh), 
                        linetype = "dashed", color = "grey50") +
    theme_publication() +
    ggplot2::labs(
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance"
    )
  
  # Add labels
  if (!is.null(label_genes)) {
    df_label <- df[df$gene %in% label_genes, ]
    p <- p + ggrepel::geom_text_repel(
      data = df_label,
      ggplot2::aes(label = gene),
      max.overlaps = 20,
      size = 3
    )
  } else if (n_top_label > 0) {
    # Label top genes by significance
    df_sig <- df[df$significance != "NS", ]
    df_sig <- df_sig[order(-df_sig$pval), ]
    df_label <- head(df_sig, n_top_label)
    
    p <- p + ggrepel::geom_text_repel(
      data = df_label,
      ggplot2::aes(label = gene),
      max.overlaps = 20,
      size = 3
    )
  }
  
  return(p)
}

#' Create paired comparison boxplot
#'
#' @param df Data frame with values
#' @param x_col Column for x-axis groups
#' @param y_col Column for y-axis values
#' @param pair_col Column for pairing
#' @param fill_col Column for fill color (optional)
#' @param log_scale Use log10 y-axis
#' @param title Plot title
#' @return ggplot object
#' @export
plot_paired_boxplot <- function(df, x_col, y_col, pair_col,
                                 fill_col = NULL, log_scale = TRUE,
                                 title = NULL) {
  
  df$x <- df[[x_col]]
  df$y <- df[[y_col]]
  df$pair <- df[[pair_col]]
  
  if (is.null(fill_col)) {
    df$fill <- df$x
  } else {
    df$fill <- df[[fill_col]]
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = fill), alpha = 0.5, 
                          outlier.shape = NA) +
    ggplot2::geom_line(ggplot2::aes(group = pair), alpha = 0.4, color = "grey40") +
    ggplot2::geom_point(size = 3, shape = 21, fill = "grey70", color = "white") +
    theme_publication() +
    ggplot2::labs(
      x = x_col,
      y = y_col,
      title = title
    )
  
  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  # Add significance if comparing 2 groups
  if (length(unique(df$x)) == 2) {
    # Perform paired test
    wide_df <- tidyr::pivot_wider(
      df[, c("pair", "x", "y")],
      names_from = x,
      values_from = y
    )
    
    groups <- unique(df$x)
    test <- wilcox.test(
      wide_df[[groups[1]]], 
      wide_df[[groups[2]]], 
      paired = TRUE
    )
    p_val <- test$p.value
    
    # Add significance annotation
    y_max <- max(df$y, na.rm = TRUE) * 1.1
    
    p <- p + ggplot2::annotate(
      "text",
      x = 1.5,
      y = y_max,
      label = sprintf("p = %.3g", p_val),
      size = 4
    )
  }
  
  return(p)
}

#' Create multi-panel figure grid
#'
#' @param plot_list List of ggplot objects
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param labels Panel labels (e.g., c("A", "B", "C"))
#' @param rel_widths Relative widths of columns
#' @param rel_heights Relative heights of rows
#' @return cowplot grid
#' @export
create_figure_panel <- function(plot_list, ncol = NULL, nrow = NULL,
                                 labels = "AUTO", 
                                 rel_widths = NULL, rel_heights = NULL) {
  
  cowplot::plot_grid(
    plotlist = plot_list,
    ncol = ncol,
    nrow = nrow,
    labels = labels,
    label_size = 14,
    label_fontface = "bold",
    rel_widths = rel_widths,
    rel_heights = rel_heights,
    align = "hv",
    axis = "tblr"
  )
}

#' Create violin plot for gene expression
#'
#' @param object Seurat object
#' @param features Genes to plot
#' @param group_by Grouping variable
#' @param split_by Split violins by this variable (optional)
#' @param pt_size Point size (0 to hide)
#' @return ggplot object
#' @export
plot_violin <- function(object, features, group_by = "seurat_clusters",
                         split_by = NULL, pt_size = 0) {
  
  # Get expression data
  df <- Seurat::FetchData(object, vars = c(features, group_by, split_by))
  
  # Reshape to long format
  df_long <- tidyr::pivot_longer(
    df,
    cols = all_of(features),
    names_to = "gene",
    values_to = "expression"
  )
  
  df_long$group <- df_long[[group_by]]
  
  p <- ggplot2::ggplot(df_long, 
                       ggplot2::aes(x = group, y = expression, fill = group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::facet_wrap(~gene, scales = "free_y") +
    theme_publication(legend = FALSE, x_axis_rotation = 45) +
    ggplot2::labs(x = group_by, y = "Expression")
  
  if (pt_size > 0) {
    p <- p + ggplot2::geom_jitter(width = 0.2, size = pt_size, alpha = 0.3)
  }
  
  if (!is.null(split_by)) {
    df_long$split <- df_long[[split_by]]
    p <- p + ggplot2::aes(fill = split)
  }
  
  return(p)
}
