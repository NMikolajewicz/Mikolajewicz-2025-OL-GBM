#!/usr/bin/env Rscript
# ==============================================================================
# main.R - GBM scRNA-seq Analysis Pipeline
# ==============================================================================
#
# Input:  Raw Seurat objects from sci-RNA-seq3 data
# Output: Annotated Seurat objects with tumor/normal classifications,
#         cell type annotations.
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Clear environment
rm(list = ls())
gc()

# Start timer
start_time <- Sys.time()

# Load configuration
source("R/utils.R")
config <- load_config("config/config.yml")

# Set seed for reproducibility
set.seed(config$reproducibility$seed)

# Setup output directories
setup_directories(config)

# Load required packages
packages <- c(
  "Seurat", "sctransform", "glmGamPoi",
  "dplyr", "tidyr", "ggplot2", "cowplot",
  "presto", "mclust", "msigdbr",
  "Matrix", "future"
)

log_message("Loading packages...")
for (pkg in packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Source analysis functions
source("R/qc.R")
source("R/integration.R")
source("R/cnv.R")
source("R/differential_expression.R")

# Configure parallel processing
options(future.globals.maxSize = config$compute$future_globals_maxsize * 1024^2)
plan("multicore", workers = config$compute$n_cores)

log_message("Setup complete.")
print_session_info(config)

# ==============================================================================
# DATA LOADING
# ==============================================================================

log_message("\n=== LOADING DATA ===")

# Define input files (update paths as needed)
# each componentn file represented separate sciRNA-seq3 experiementaln run
input_files <- list(
  p9 = file.path(config$paths$raw_data_dir, "p9_hs_pr_gbm.rds"),
  p14 = file.path(config$paths$raw_data_dir, "p14_hs_pr_gbm.rds"),
  p16 = file.path(config$paths$raw_data_dir, "p16_hs_pr_gbm.rds")
)

# Check files exist
for (name in names(input_files)) {
  if (!file.exists(input_files[[name]])) {
    log_message(sprintf("WARNING: Input file not found: %s", 
                        input_files[[name]]), "WARNING")
  }
}

# Load Seurat objects (adjust based on actual file format)
so_list <- list()

for (name in names(input_files)) {
  if (file.exists(input_files[[name]])) {
    log_message(sprintf("Loading %s...", name))
    
    if (grepl("\\.rds$", input_files[[name]], ignore.case = TRUE)) {
      so_list[[name]] <- readRDS(input_files[[name]])
    } else if (grepl("\\.Rdata$", input_files[[name]], ignore.case = TRUE)) {
      env <- new.env()
      load(input_files[[name]], envir = env)
      so_list[[name]] <- get("so", envir = env)
    }
    
    log_message(sprintf("  Loaded: %d cells, %d genes", 
                        ncol(so_list[[name]]), nrow(so_list[[name]])))
  }
}

# ==============================================================================
# PREPROCESSING
# ==============================================================================

log_message("\n=== PREPROCESSING ===")

# Process each sample
so_list_processed <- list()

for (name in names(so_list)) {
  log_message(sprintf("\nProcessing sample: %s", name))
  
  # Preprocess (QC, normalize, reduce dimensions)
  so_list_processed[[name]] <- preprocess_sample(
    so_list[[name]], 
    config, 
    sample_name = name
  )
  
  # Add sample metadata
  so_list_processed[[name]]$batch <- name
  so_list_processed[[name]]$sample <- name
}

# ==============================================================================
# INTEGRATION
# ==============================================================================

log_message("\n=== INTEGRATION ===")

# Integrate samples using rPCA (as per manuscript)
so_integrated <- integrate_rpca(so_list_processed, config)

# Rename samples using manuscript labels
sample_map <- config$sample_metadata$rename_barcode
idh_map <- config$sample_metadata$idh_status

so_integrated$sample_id <- sample_map[so_integrated$batch]
so_integrated$idh_status <- idh_map[so_integrated$sample_id]

log_message(sprintf("Integrated object: %d cells", ncol(so_integrated)))

# Apply BBKNN for additional batch correction visualization
so_integrated <- integrate_bbknn(so_integrated, batch_key = "sample_id", config)

# ==============================================================================
# TUMOR/NORMAL CLASSIFICATION
# ==============================================================================

log_message("\n=== TUMOR CALLING ===")

# Identify microglia clusters for CNV reference
# clusters 3, 5, 16, 13, 15
microglia_clusters <- config$cnv$reference_clusters

# Run CNV-based tumor calling
so_integrated <- call_tumor_cells(
  so_integrated, 
  config, 
  microglia_clusters = microglia_clusters
)

# Summarize results
tumor_summary <- data.frame(
  sample = unique(so_integrated$sample_id),
  n_tumor = tapply(so_integrated$is_tumor, 
                   so_integrated$sample_id, sum, na.rm = TRUE),
  n_normal = tapply(!so_integrated$is_tumor, 
                    so_integrated$sample_id, sum, na.rm = TRUE)
)
tumor_summary$pct_tumor <- tumor_summary$n_tumor / 
  (tumor_summary$n_tumor + tumor_summary$n_normal) * 100

log_message("Tumor calling summary:")
print(tumor_summary)

# Save tumor summary table
save_table(tumor_summary, "TableS_tumor_classification_summary", config)

# ==============================================================================
# CELL TYPE ANNOTATION
# ==============================================================================

log_message("\n=== CELL TYPE ANNOTATION ===")

# Score cell type markers
cell_type_markers <- config$cell_type_markers

for (ct_name in names(cell_type_markers)) {
  markers <- cell_type_markers[[ct_name]]
  markers_present <- intersect(markers, rownames(so_integrated))
  
  if (length(markers_present) >= 2) {
    so_integrated <- Seurat::AddModuleScore(
      so_integrated,
      features = list(markers_present),
      name = paste0(ct_name, "_score"),
      ctrl = 100
    )
    log_message(sprintf("  Scored %s: %d/%d markers", 
                        ct_name, length(markers_present), length(markers)))
  }
}

# Assign cell types based on highest score
# Focus on non-tumor cells
score_cols <- grep("_score1$", colnames(so_integrated@meta.data), value = TRUE)
cell_type_names <- gsub("_score1$", "", score_cols)

# Get max score per cell
score_mat <- so_integrated@meta.data[, score_cols]
max_score_idx <- apply(score_mat, 1, which.max)
so_integrated$cell_type_pred <- cell_type_names[max_score_idx]

# Set tumor cells
so_integrated$cell_type_final <- ifelse(
  so_integrated$is_tumor, 
  "tumor", 
  so_integrated$cell_type_pred
)

# Summary
ct_summary <- data.frame(table(so_integrated$cell_type_final))
colnames(ct_summary) <- c("cell_type", "n_cells")
ct_summary$pct <- ct_summary$n_cells / sum(ct_summary$n_cells) * 100

log_message("Cell type summary:")
print(ct_summary)

# Save cell type summary
save_table(ct_summary, "TableS_celltype_summary", config)

# ==============================================================================
# SUBSET ANALYSIS
# ==============================================================================

log_message("\n=== SUBSET ANALYSIS ===")

# Split into tumor and normal populations
so_tumor <- subset(so_integrated, subset = is_tumor == TRUE)
so_normal <- subset(so_integrated, subset = is_tumor == FALSE)

log_message(sprintf("Tumor cells: %d", ncol(so_tumor)))
log_message(sprintf("Normal cells: %d", ncol(so_normal)))

# Re-embed tumor cells
so_tumor <- re_embed_subset(so_tumor, config, re_cluster = TRUE)

# Re-embed normal cells
so_normal <- re_embed_subset(so_normal, config, re_cluster = TRUE)

# ==============================================================================
# DIFFERENTIAL EXPRESSION
# ==============================================================================

log_message("\n=== DIFFERENTIAL EXPRESSION ===")

# DE for normal cell types
de_normal <- run_de_presto(so_normal, group_by = "cell_type_final", config)
top_markers_normal <- get_top_markers(de_normal, n_top = 20)

save_table(top_markers_normal, "TableS_normal_celltype_markers", config)

# DE for tumor clusters
de_tumor <- run_de_presto(so_tumor, group_by = "seurat_clusters", config)
top_markers_tumor <- get_top_markers(de_tumor, n_top = 20)

save_table(top_markers_tumor, "TableS_tumor_cluster_markers", config)

# ==============================================================================
# ABUNDANCE ANALYSIS
# ==============================================================================

log_message("\n=== ABUNDANCE ANALYSIS ===")

# Calculate cell type abundances
df_abundance <- calculate_abundance(
  so_integrated, 
  cell_type_col = "cell_type_final",
  sample_col = "sample_id"
)

# Add primary/recurrent info
df_abundance$PR <- ifelse(grepl("^P", df_abundance$sample), "Primary", "Recurrent")
df_abundance$pair_id <- gsub("^P|^R", "", df_abundance$sample)

# Compare primary vs recurrent
abundance_comparison <- compare_abundance(
  df_abundance,
  condition_col = "PR",
  conditions = c("Primary", "Recurrent"),
  pair_col = "pair_id"
)

log_message("Abundance comparison (Primary vs Recurrent):")
print(abundance_comparison[order(abundance_comparison$padj_paired), ])

save_table(abundance_comparison, "TableS_abundance_PR_comparison", config)
save_table(df_abundance, "TableS_celltype_abundances", config)

# ==============================================================================
# GENERATE FIGURES
# ==============================================================================

log_message("\n=== GENERATING FIGURES ===")

# Figure 1A: Cell type fractional abundance
plot_list <- list()

for (ct in unique(df_abundance$cell_type)) {
  plot_list[[ct]] <- plot_abundance_comparison(
    df_abundance, ct, "PR", "pair_id"
  )
}

fig1a <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)
save_figure(fig1a, "Fig1A_fractional_abundance.png", 
            width = 12, height = 8, config = config)

# Figure S1: UMAP by cell type and sample
fig_s1a <- Seurat::DimPlot(so_integrated, reduction = "umap", 
                           group.by = "cell_type_final", 
                           label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("Cell Types")

fig_s1b <- Seurat::DimPlot(so_integrated, reduction = "umap",
                           group.by = "sample_id") +
  ggplot2::ggtitle("Samples")

fig_s1c <- Seurat::DimPlot(so_integrated, reduction = "umap",
                           group.by = "is_tumor") +
  ggplot2::ggtitle("Tumor Status")

fig_s1 <- cowplot::plot_grid(fig_s1a, fig_s1b, fig_s1c, ncol = 3)
save_figure(fig_s1, "FigS1_UMAP_overview.png", 
            width = 18, height = 5, config = config)

# Cell type marker dot plot
marker_genes <- c(
  # Endothelial
  "CD34", "CDH5",
  # Astrocyte
  "AGT", "GFAP",
  # Oligodendrocyte
  "PLP1", "MOG",
  # OPC
  "PDGFRA", "PTPRZ1",
  # Microglia
  "P2RY12", "CD163",
  # Neurons
  "RBFOX3", "CELF4",
  # T cells
  "CD3D", "CD3E",
  # Tumor markers
  "EGFR", "SOX2", "PROM1"
)

marker_genes_present <- intersect(marker_genes, rownames(so_normal))

df_dot_normal <- prepare_dotplot_data(
  so_normal, 
  features = marker_genes_present,
  group_by = "cell_type_final"
)

fig_dot <- plot_markers_dot(df_dot_normal, feature_order = marker_genes_present)
save_figure(fig_dot, "FigS1_marker_dotplot.png",
            width = 8, height = 10, config = config)

# ==============================================================================
# SAVE PROCESSED DATA
# ==============================================================================

log_message("\n=== SAVING DATA ===")

# Save integrated object
saveRDS(so_integrated, 
        file.path(config$paths$processed_data_dir, "so_integrated.rds"))

# Save tumor and normal subsets
saveRDS(so_tumor,
        file.path(config$paths$processed_data_dir, "so_tumor.rds"))
saveRDS(so_normal,
        file.path(config$paths$processed_data_dir, "so_normal.rds"))

# Final session info
sessionInfo()
