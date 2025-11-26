#===============================================================================
# Cytokine Profiling and Ligand-Receptor Prevalence Analysis
#===============================================================================
# 
# Description: Consolidated script for cytokine analysis and ligand-receptor
#              interaction prevalence calculations for the OL-GBM manuscript.
#
# Purpose: cytokine profiling and LR prevalence metrics
#          from scRNA-seq data across GBM cohorts.
# 
# Author: N. Mikolajewicz
#
# R version used: R >= 4.1.0
# Key package versions:
#   - Seurat >= 4.0.0
#   - dplyr >= 1.0.0
#   - tidyr >= 1.0.0
#   - readxl >= 1.4.0
#   - circlize >= 0.4.16
#
#===============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Clear environment
rm(list = ls())
invisible(gc())

# Set seed for reproducibility
set.seed(12345)

# Record start time
start_time <- proc.time()

# List of required packages
packages_required <- c(
  "Seurat",
  "sctransform",
  "dplyr",
  "tidyr",
  "readxl",
  "ggplot2",
  "ggrepel",
  "RColorBrewer",
  "reshape2",
  "circlize",
  "cowplot"
)

# Load packages (install if necessary)
for (pkg in packages_required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing package: ", pkg))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Load scMiko package (custom functions - must be installed separately)
# This package contains helper functions including:
#   - theme_miko(), scale_color_miko(), scale_fill_miko()
#   - p2z(), z2p() - p-value to z-score conversions
#   - LR.db - ligand-receptor database
#   - categoricalColPal() - color palette generator
if (!requireNamespace("scMiko", quietly = TRUE)) {
  stop("The scMiko package is required. Install from: github.com/Mikolajewicz-2025-OL-GBM")
}
library(scMiko)

# Set global options
options(stringsAsFactors = FALSE)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Two-sample Z-test
#' @description Compares means of two samples using z-test approximation
#' @param x First sample
#' @param y Second sample
#' @param return.p Return p-value (TRUE) or z-statistic (FALSE)
#' @param winsorize Winsorization threshold for extreme values
#' @return P-value or z-statistic
z.test.2s <- function(x, y, return.p = TRUE, winsorize = 30) {
  x.mean <- mean(x, na.rm = TRUE)
  y.mean <- mean(y, na.rm = TRUE)
  
  if (length(x) > 1) {
    x.sd <- sd(x, na.rm = TRUE)
    x.se <- x.sd / sqrt(length(x))
  } else {
    x.sd <- x.se <- 0
  }
  
  if (length(y) > 1) {
    y.sd <- sd(y, na.rm = TRUE)
    y.se <- y.sd / sqrt(length(y))
  } else {
    y.sd <- y.se <- 0
  }
  
  delta <- y.mean - x.mean
  se <- sqrt((x.se^2) + (y.se^2))
  z <- delta / se
  
  # Winsorize extreme values
  tryCatch({
    if (z > winsorize) z <- winsorize
    if (z < (-1 * winsorize)) z <- -winsorize
  }, error = function(e) {})
  
  if (return.p) {
    return(z2p(z))
  } else {
    return(z)
  }
}

#' One-sample Z-test
#' @description Tests if sample mean differs from zero
#' @param x Sample
#' @param return.p Return p-value (TRUE) or z-statistic (FALSE)
#' @param winsorize Winsorization threshold for extreme values
#' @return P-value or z-statistic
z.test.1s <- function(x, return.p = TRUE, winsorize = 30) {
  x.mean <- mean(x, na.rm = TRUE)
  
  if (length(x) > 1) {
    x.sd <- sd(x, na.rm = TRUE)
    x.se <- x.sd / sqrt(length(x))
    z <- x.mean / x.se
  } else {
    x.sd <- x.se <- 0
    z <- NA
  }
  
  if (!is.na(z)) {
    if (z > winsorize) z <- winsorize
    if (z < (-1 * winsorize)) z <- -winsorize
    p <- z2p(z)
  } else {
    p <- NA
  }
  
  if (return.p) {
    return(p)
  } else {
    return(z)
  }
}

#' Extract LR expression data from averaged expression matrices
#' @description Helper function to extract ligand and receptor expression 
#'              for a given cohort and condition
#' @param dat Averaged expression matrix
#' @param cohort Cohort name
#' @param PR Primary ("P") or Recurrent ("R") status
#' @param source Source cell type (ligand-expressing)
#' @param target Target cell type (receptor-expressing)
#' @param LIG Ligand gene symbol (from parent environment)
#' @param REC Receptor gene symbol (from parent environment)
#' @return Data frame with LR expression data
getData <- function(dat, cohort, PR, source = "oligodendrocyte", target = "tumor") {
  data.names <- colnames(dat)
  data.names <- unique(gsub("_", "", gsub(paste0(source, "|", target), "", data.names)))
  
  dat.new.all <- NULL
  for (i in seq_along(data.names)) {
    source.name <- paste0(source, "_", data.names[i])
    target.name <- paste0(target, "_", data.names[i])
    
    if ((LIG %in% rownames(dat) & REC %in% rownames(dat)) & 
        (source.name %in% colnames(dat)) & 
        (target.name %in% colnames(dat))) {
      LIG_EXPR <- as.numeric(dat[LIG, source.name])
      REC_EXPR <- as.numeric(dat[REC, target.name])
      dat.new <- data.frame(
        cohort = cohort,
        PR = PR,
        patient = data.names[i],
        LIG = LIG,
        REC = REC,
        LIG_EXPR = LIG_EXPR,
        REC_EXPR = REC_EXPR
      )
    } else {
      dat.new <- NULL
    }
    dat.new.all <- bind_rows(dat.new.all, dat.new)
  }
  return(dat.new.all)
}

# ==============================================================================
# CONFIGURATION - FILE PATHS
# ==============================================================================

# Input files - modify these paths according to your directory structure
# All paths should be relative to the project root

config <- list(
  # Cytokine profiling data (96-plex assay results)
  cytokine_data_path = "data/cytokine_analysis_280224.xlsx",
  
  # Preprocessed Seurat objects containing tumor and oligodendrocyte cells
  # from Mikolajewicz, Abdelfattah, and Wang cohorts
  seurat_data_path = "data/PR4_preprocessed_tumor_and_other_miko_wang_abdel_140623.rds",
  
  # Output directory
  output_dir = "outputs/PR14_results/"
)

# Create output directory if it doesn't exist
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
}

# ==============================================================================
# SECTION 1: LOAD AND PREPROCESS CYTOKINE DATA
# ==============================================================================

message("=== SECTION 1: Loading cytokine profiling data ===")

# Load cytokine profiling data from 96-plex assay
df.dat <- readxl::read_xlsx(
  path = config$cytokine_data_path,
  sheet = "clean_data"
)

# Extract cytokine names (all columns except metadata)
cytokines.all <- colnames(df.dat)
cytokines.all <- cytokines.all[!(cytokines.all %in% c("ID", "Sample1", "Sample2", "Sample3"))]

# Clean data matrix - replace "OOR <" (out of range) values with 0
df.mat <- as.matrix(df.dat[, cytokines.all])
df.mat[df.mat %in% "OOR <"] <- 0
df.mat[df.mat %in% "OOR"] <- 0
df.dat[, cytokines.all] <- df.mat

# Convert to long format for analysis
df.dat.long <- pivot_longer(data = df.dat, cols = all_of(cytokines.all))
colnames(df.dat.long) <- c("ID", "Sample1", "Sample2", "Sample3", "cytokine", "expr")
df.dat.long$expr <- as.numeric(df.dat.long$expr)

message(paste0("  Loaded ", length(cytokines.all), " cytokines from ", 
               nrow(df.dat), " samples"))

# ==============================================================================
# SECTION 2: OL CYTOKINE ANALYSIS
# ==============================================================================

message("=== SECTION 2: OPC cytokine analysis ===")

# Filter for OPC and NCC (negative control cell) samples
df.opc <- df.dat.long %>% dplyr::filter(Sample1 %in% c("OPC", "NCC"))

# Calculate summary statistics comparing OPC to NCC
# Using z-test to assess significance of OPC cytokine expression
# relative to NCC baseline
df.opc.sum <- df.opc %>%
  dplyr::group_by(cytokine) %>%
  dplyr::summarise(
    mean.opc = mean(expr[Sample1 %in% "OPC"]),
    sd.opc = sd(expr[Sample1 %in% "OPC"]),
    mean.ncc = mean(expr[Sample1 %in% "NCC"]),
    delta = mean(expr[Sample1 %in% "OPC"]) - mean(expr[Sample1 %in% "NCC"]),
    pval = z.test.1s(expr[Sample1 %in% "OPC"] - mean(expr[Sample1 %in% "NCC"])),
    z = z.test.1s(expr[Sample1 %in% "OPC"] - mean(expr[Sample1 %in% "NCC"]), 
                  return.p = FALSE),
    .groups = "drop"
  )

# Calculate coefficient of variation
df.opc.sum$cv <- df.opc.sum$sd.opc / df.opc.sum$mean.opc

# Handle NA p-values
df.opc.sum$pval[is.na(df.opc.sum$pval)] <- 1

# Multiple testing correction (Benjamini-Hochberg)
df.opc.sum$padj <- p.adjust(df.opc.sum$pval, method = "BH")

# Recalculate z-scores from corrected p-values
df.opc.sum$z <- p2z(df.opc.sum$pval) * sign(df.opc.sum$delta)
df.opc.sum$crank <- rank(df.opc.sum$z, ties.method = "random")
df.opc.sum$logp <- -log10(df.opc.sum$padj) * sign(df.opc.sum$delta)
df.opc.sum$logp[df.opc.sum$logp < 0] <- 0

# Identify significantly upregulated cytokines in OPC
fdr.sig.opc <- 0.05
df.opc.sum.top <- df.opc.sum %>% 
  dplyr::filter(padj < fdr.sig.opc, delta > 0)

message(paste0("  Found ", nrow(df.opc.sum.top), 
               " significantly upregulated cytokines in OPC (FDR < ", 
               fdr.sig.opc, ")"))

# ==============================================================================
# SECTION 3: GBM CYTOKINE ANALYSIS
# ==============================================================================

message("=== SECTION 3: GBM cytokine analysis ===")

# Filter for GBM (pGBM + rGBM) and NCC samples
df.g <- df.dat.long %>% dplyr::filter(Sample2 %in% c("pGBM", "rGBM", "NCC"))
df.g$Sample4 <- ifelse(df.g$Sample1 %in% "NCC", "NCC", "GBM")

# Calculate GBM vs NCC summary statistics
df.g.sum <- df.g %>%
  dplyr::group_by(cytokine) %>%
  dplyr::summarise(
    mean.gbm = mean(expr[Sample4 %in% "GBM"]),
    sd.gbm = sd(expr[Sample4 %in% "GBM"]),
    mean.ncc = mean(expr[Sample4 %in% "NCC"]),
    delta = mean(expr[Sample4 %in% "GBM"]) - mean(expr[Sample4 %in% "NCC"]),
    pval = z.test.1s(expr[Sample4 %in% "GBM"] - mean(expr[Sample4 %in% "NCC"])),
    z = z.test.1s(expr[Sample4 %in% "GBM"] - mean(expr[Sample4 %in% "NCC"]), 
                  return.p = FALSE),
    .groups = "drop"
  )

# Handle NA p-values and multiple testing correction
df.g.sum$pval[is.na(df.g.sum$pval)] <- 1
df.g.sum$padj <- p.adjust(df.g.sum$pval, method = "BH")
df.g.sum$z <- p2z(df.g.sum$pval) * sign(df.g.sum$delta)
df.g.sum$crank <- rank(df.g.sum$z, ties.method = "random")
df.g.sum$logp <- -log10(df.g.sum$padj) * sign(df.g.sum$delta)
df.g.sum$logp[df.g.sum$logp < 0] <- 0

# Identify significantly expressed GBM cytokines
fdr.sig.gbm <- 0.1
df.g.sum.top <- df.g.sum %>% 
  dplyr::filter(padj < fdr.sig.gbm, delta > 0)

message(paste0("  Found ", nrow(df.g.sum.top), 
               " significantly upregulated cytokines in GBM (FDR < ", 
               fdr.sig.gbm, ")"))

# ==============================================================================
# SECTION 4: PRIMARY VS RECURRENT GBM COMPARISON
# ==============================================================================

message("=== SECTION 4: Primary vs Recurrent GBM comparison ===")

# Filter for pGBM and rGBM samples
df.gbm <- df.dat.long %>% dplyr::filter(Sample2 %in% c("pGBM", "rGBM", "NCC"))

# Calculate primary vs recurrent GBM comparison
df.gbm.sum <- df.gbm %>%
  dplyr::group_by(cytokine) %>%
  dplyr::summarise(
    mean.pGBM = mean(expr[Sample2 %in% "pGBM"]) - mean(expr[Sample1 %in% "NCC"]),
    sd.pGBM = sd(expr[Sample2 %in% "pGBM"] - mean(expr[Sample1 %in% "NCC"])),
    mean.rGBM = mean(expr[Sample2 %in% "rGBM"]) - mean(expr[Sample1 %in% "NCC"]),
    sd.rGBM = sd(expr[Sample2 %in% "rGBM"] - mean(expr[Sample1 %in% "NCC"])),
    pval1 = z.test.2s(
      expr[Sample2 %in% "pGBM"] - mean(expr[Sample1 %in% "NCC"]),
      expr[Sample2 %in% "rGBM"] - mean(expr[Sample1 %in% "NCC"]),
      return.p = TRUE, winsorize = 5
    ),
    pval2 = z.test.1s(
      log10(expr[Sample2 %in% "pGBM"] + 1) - log10(expr[Sample2 %in% "rGBM"] + 1),
      return.p = TRUE, winsorize = 5
    ),
    .groups = "drop"
  )

# Handle NA values and select minimum p-value
df.gbm.sum$pval1[is.na(df.gbm.sum$pval1)] <- 1
df.gbm.sum$pval2[is.na(df.gbm.sum$pval2)] <- 1
df.gbm.sum$pval <- df.gbm.sum$pval1
df.gbm.sum$pval[df.gbm.sum$pval2 < df.gbm.sum$pval1] <- 
  df.gbm.sum$pval2[df.gbm.sum$pval2 < df.gbm.sum$pval1]

# Ensure non-negative means
df.gbm.sum$mean.pGBM[df.gbm.sum$mean.pGBM < 0] <- 0
df.gbm.sum$mean.rGBM[df.gbm.sum$mean.rGBM < 0] <- 0

# Calculate delta and log fold change
df.gbm.sum$delta <- df.gbm.sum$mean.rGBM - df.gbm.sum$mean.pGBM
df.gbm.sum$standardized_dif <- df.gbm.sum$delta / 
  sqrt((df.gbm.sum$sd.pGBM^2) + df.gbm.sum$sd.rGBM^2)
df.gbm.sum$lfc <- log10(df.gbm.sum$mean.rGBM + 1) - 
  log10(df.gbm.sum$mean.pGBM + 1)
df.gbm.sum$padj <- p.adjust(df.gbm.sum$pval, method = "BH")
df.gbm.sum$z <- p2z(df.gbm.sum$pval) * sign(df.gbm.sum$delta)

# Identify significantly differentially expressed cytokines
sig.thresh.pr <- 0.2
df.gbm.sum.sig <- df.gbm.sum %>%
  dplyr::filter(cytokine %in% df.g.sum.top$cytokine | padj < sig.thresh.pr)
df.gbm.sum.top <- df.gbm.sum.sig %>% 
  dplyr::filter(cytokine %in% df.g.sum.top$cytokine | pval < sig.thresh.pr)

message(paste0("  Found ", nrow(df.gbm.sum.top), 
               " differentially expressed cytokines between pGBM and rGBM"))

# Identify which GBM cytokines are actually expressed
df.gbm.sum$z.pGBM <- df.gbm.sum$mean.pGBM / (df.gbm.sum$sd.pGBM / sqrt(3))
df.gbm.sum$z.rGBM <- df.gbm.sum$mean.rGBM / (df.gbm.sum$sd.rGBM / sqrt(3))
df.gbm.sum$q.pGBM <- p.adjust(z2p(df.gbm.sum$z.pGBM), method = "BH")
df.gbm.sum$q.rGBM <- p.adjust(z2p(df.gbm.sum$z.rGBM), method = "BH")

qthresh <- 0.2
which.gbm.expressed <- (df.gbm.sum %>%
                          dplyr::filter(q.pGBM < qthresh | q.rGBM < qthresh))$cytokine

# ==============================================================================
# SECTION 5: LOAD LIGAND-RECEPTOR DATABASE AND CYTOKINE MAPPING
# ==============================================================================

message("=== SECTION 5: Loading LR database and cytokine mapping ===")

# Load cytokine to gene name mapping
df.ligs <- readxl::read_xlsx(
  path = config$cytokine_data_path,
  sheet = "ligands"
)

# Create protein-to-gene mapping vector
prot2gene <- df.ligs$Ligand_Gene
names(prot2gene) <- df.ligs$Cytokine

# Load LR database from scMiko package
# This database combines: SCSR, iTALK, Fantom5, NATMI
data("LR.db")
df.LR <- LR.db[["LR.df"]]

# Subset LR database to cytokines profiled in the 96-plex assay
df.LR.subset <- df.LR %>% 
  dplyr::filter(Ligand %in% df.ligs$Ligand_Gene)

message(paste0("  Loaded ", nrow(df.LR.subset), 
               " LR pairs for profiled cytokines"))

# ==============================================================================
# SECTION 6: LOAD AND PREPROCESS scRNA-seq DATA
# ==============================================================================

message("=== SECTION 6: Loading scRNA-seq data ===")

# Load preprocessed Seurat objects for three GBM cohorts
so.all <- readRDS(config$seurat_data_path)

# Subset to tumor and oligodendrocyte cells only
so.miko <- so.all[["mikolajewicz"]][
  , so.all[["mikolajewicz"]]@meta.data[["class_consensus"]] %in% 
    c("tumor", "oligodendrocyte")
]
so.abdel <- so.all[["abdelfattah"]][
  , so.all[["abdelfattah"]]@meta.data[["class_consensus"]] %in% 
    c("tumor", "oligodendrocyte")
]
so.wang <- so.all[["wang"]][
  , so.all[["wang"]]@meta.data[["class_consensus"]] %in% 
    c("tumor", "oligodendrocyte")
]

# Clear memory
rm(so.all)
invisible(gc())

message(paste0("  Mikolajewicz cohort: ", ncol(so.miko), " cells"))
message(paste0("  Abdelfattah cohort: ", ncol(so.abdel), " cells"))
message(paste0("  Wang cohort: ", ncol(so.wang), " cells"))

# ==============================================================================
# SECTION 7: CALCULATE AVERAGE EXPRESSION FOR LR ANALYSIS
# ==============================================================================

message("=== SECTION 7: Calculating average expression per cell type ===")

# Define assay to use
which.assay <- "SCT"
LR.features <- unique(c(df.LR.subset$Ligand, df.LR.subset$Receptor))

# --- Mikolajewicz cohort ---
df.miko.expr <- AverageExpression(
  object = so.miko,
  assays = which.assay,
  features = LR.features,
  return.seurat = FALSE,
  group.by = c("class_consensus", "sample"),
  slot = "data",
  verbose = FALSE
)[[1]]
df.miko.expr <- as.data.frame(df.miko.expr)
df.miko.expr.p <- df.miko.expr[, grepl("_P", colnames(df.miko.expr))]
df.miko.expr.r <- df.miko.expr[, grepl("_R", colnames(df.miko.expr))]

# --- Abdelfattah cohort ---
# Create unique patient identifiers
so.abdel@meta.data$PR2 <- so.abdel@meta.data$PR
usample <- unique(so.abdel@meta.data$sample)
usample2 <- seq_along(usample)
names(usample2) <- usample
so.abdel@meta.data$PR2 <- paste0(
  so.abdel@meta.data$PR2, 
  usample2[so.abdel@meta.data$sample]
)

df.abdel.expr <- AverageExpression(
  object = so.abdel,
  assays = which.assay,
  features = LR.features,
  return.seurat = FALSE,
  group.by = c("class_consensus", "PR2"),
  slot = "data",
  verbose = FALSE
)[[1]]
df.abdel.expr <- as.data.frame(df.abdel.expr)
df.abdel.expr.p <- df.abdel.expr[, grepl("_P", colnames(df.abdel.expr))]
df.abdel.expr.r <- df.abdel.expr[, grepl("_R", colnames(df.abdel.expr))]

# --- Wang cohort ---
so.wang@meta.data$PR2 <- so.wang@meta.data$PR
so.wang@meta.data$PR2[grepl("Primary", so.wang@meta.data$PR)] <- "P"
so.wang@meta.data$PR2[grepl("Recurrent", so.wang@meta.data$PR)] <- "R"
so.wang@meta.data$PR2 <- paste0(so.wang@meta.data$PR2, so.wang@meta.data$pair)

df.wang.expr <- AverageExpression(
  object = so.wang,
  assays = which.assay,
  features = LR.features,
  return.seurat = FALSE,
  group.by = c("class_consensus", "PR2"),
  slot = "data",
  verbose = FALSE
)[[1]]
df.wang.expr <- as.data.frame(df.wang.expr)
df.wang.expr.p <- df.wang.expr[, grepl("_P", colnames(df.wang.expr))]
df.wang.expr.r <- df.wang.expr[, grepl("_R", colnames(df.wang.expr))]

message("  Calculated average expression matrices for all cohorts")

# ==============================================================================
# SECTION 8: CALCULATE LR PREVALENCE - OLIGODENDROCYTE SOURCE
# ==============================================================================

message("=== SECTION 8: Calculating LR prevalence (OL source -> GBM target) ===")

# Map OPC-expressed cytokines to gene symbols
df.opc.sum.top$LIG <- prot2gene[df.opc.sum.top$cytokine]

# Filter LR database to OPC-expressed ligands
df.LR.olig <- df.LR.subset %>% 
  dplyr::filter(Ligand %in% df.opc.sum.top$LIG)

# Calculate LR interactions for all pairs across all cohorts
df.olig.interactions <- NULL

for (i in seq_len(nrow(df.LR.olig))) {
  LIG <- df.LR.olig$Ligand[i]
  REC <- df.LR.olig$Receptor[i]
  
  # Collect interactions from all cohorts and conditions
  df.olig.interactions <- bind_rows(
    df.olig.interactions,
    getData(df.miko.expr.p, "Mikolajewicz", "P"),
    getData(df.miko.expr.r, "Mikolajewicz", "R"),
    getData(df.abdel.expr.p, "Abdelfattah", "P"),
    getData(df.abdel.expr.r, "Abdelfattah", "R"),
    getData(df.wang.expr.p, "Wang", "P"),
    getData(df.wang.expr.r, "Wang", "R")
  )
}

# Calculate LR interaction score (ligand × receptor expression product)
df.olig.interactions$LR <- df.olig.interactions$LIG_EXPR * 
  df.olig.interactions$REC_EXPR
df.olig.interactions$LR_name <- paste0(
  df.olig.interactions$LIG, "_", 
  df.olig.interactions$REC
)

# Rank LR interactions within each patient
df.olig.interactions <- df.olig.interactions %>%
  dplyr::group_by(cohort, PR, patient) %>%
  dplyr::mutate(LR_rank = rank(LR, ties.method = "random"))

# Calculate LR prevalence (P_LR): fraction of samples with active interaction
# An interaction is "active" if L × R > 0
df.olig.interactions.sum <- df.olig.interactions %>%
  dplyr::group_by(LR_name, LIG, REC) %>%
  dplyr::summarise(
    # Primary GBM statistics
    LR.mean.P = mean(LR[PR %in% "P"], na.rm = TRUE),
    LR.prop.p = mean(LR[PR %in% "P"] > 0, na.rm = TRUE),  # P_LR for primary
    LIG.EXPR.P = mean(LIG_EXPR[PR %in% "P"], na.rm = TRUE),
    REC.EXPR.P = mean(REC_EXPR[PR %in% "P"], na.rm = TRUE),
    
    # Recurrent GBM statistics
    LR.mean.R = mean(LR[PR %in% "R"], na.rm = TRUE),
    LR.prop.R = mean(LR[PR %in% "R"] > 0, na.rm = TRUE),  # P_LR for recurrent
    LIG.EXPR.R = mean(LIG_EXPR[PR %in% "R"], na.rm = TRUE),
    REC.EXPR.R = mean(REC_EXPR[PR %in% "R"], na.rm = TRUE),
    
    # Combined statistics
    delta = mean(LR[PR %in% "R"], na.rm = TRUE) - mean(LR[PR %in% "P"], na.rm = TRUE),
    LR.prop.p.all = mean(LR[PR %in% c("P", "R")] > 0, na.rm = TRUE),  # Overall P_LR
    
    # Sample counts
    n.p.active = sum(LR[PR %in% "P"] > 0, na.rm = TRUE),
    n.p.total = sum(PR %in% "P", na.rm = TRUE),
    n.r.active = sum(LR[PR %in% "R"] > 0, na.rm = TRUE),
    n.r.total = sum(PR %in% "R", na.rm = TRUE),
    
    # Statistical tests
    delta.prop.pval = suppressWarnings(
      prop.test(x = c(sum(LR[PR %in% "P"] > 0), sum(LR[PR %in% "R"] > 0)),
                n = c(sum(PR %in% "P"), sum(PR %in% "R")))[["p.value"]]
    ),
    delta.pval = suppressWarnings(
      wilcox.test(LR[PR %in% "P"], LR[PR %in% "R"])$p.value
    ),
    P.pval = suppressWarnings(wilcox.test(LR[PR %in% "P"])$p.value),
    R.pval = suppressWarnings(wilcox.test(LR[PR %in% "R"])$p.value),
    LFC = log(mean(LR[PR %in% "R"]) + 1) - log(mean(LR[PR %in% "P"]) + 1),
    LR.rank.P = mean(LR_rank[PR %in% "P"]),
    LR.rank.R = mean(LR_rank[PR %in% "R"]),
    LR.rank.pval = suppressWarnings(
      wilcox.test(LR_rank[PR %in% "P"], LR_rank[PR %in% "R"])$p.value
    ),
    .groups = "drop"
  )

# Calculate additional metrics and multiple testing correction
df.olig.interactions.sum$prop.delta <- df.olig.interactions.sum$LR.prop.R - 
  df.olig.interactions.sum$LR.prop.p
df.olig.interactions.sum$delta.fdr <- p.adjust(
  df.olig.interactions.sum$delta.pval, method = "BH"
)
df.olig.interactions.sum$delta.prop.fdr <- p.adjust(
  df.olig.interactions.sum$delta.prop.pval, method = "BH"
)
df.olig.interactions.sum$P.padj <- p.adjust(
  df.olig.interactions.sum$P.pval, method = "BH"
)
df.olig.interactions.sum$R.padj <- p.adjust(
  df.olig.interactions.sum$R.pval, method = "BH"
)

message(paste0("  Calculated prevalence for ", 
               nrow(df.olig.interactions.sum), " OL-source LR pairs"))

# ==============================================================================
# SECTION 9: CALCULATE LR PREVALENCE - GBM SOURCE
# ==============================================================================

message("=== SECTION 9: Calculating LR prevalence (GBM source -> OL target) ===")

# Map GBM-expressed cytokines to gene symbols
df.gbm.sum.top$LIG <- prot2gene[df.gbm.sum.top$cytokine]

# Filter LR database to GBM-expressed ligands
df.LR.tumor <- df.LR.subset %>% 
  dplyr::filter(Ligand %in% df.gbm.sum.top$LIG)

# Calculate LR interactions for all pairs across all cohorts
df.tumor.interactions <- NULL

for (i in seq_len(nrow(df.LR.tumor))) {
  LIG <- df.LR.tumor$Ligand[i]
  REC <- df.LR.tumor$Receptor[i]
  
  # Collect interactions (source = tumor, target = oligodendrocyte)
  df.tumor.interactions <- bind_rows(
    df.tumor.interactions,
    getData(df.miko.expr.p, "Mikolajewicz", "P", 
            source = "tumor", target = "oligodendrocyte"),
    getData(df.miko.expr.r, "Mikolajewicz", "R", 
            source = "tumor", target = "oligodendrocyte"),
    getData(df.abdel.expr.p, "Abdelfattah", "P", 
            source = "tumor", target = "oligodendrocyte"),
    getData(df.abdel.expr.r, "Abdelfattah", "R", 
            source = "tumor", target = "oligodendrocyte"),
    getData(df.wang.expr.p, "Wang", "P", 
            source = "tumor", target = "oligodendrocyte"),
    getData(df.wang.expr.r, "Wang", "R", 
            source = "tumor", target = "oligodendrocyte")
  )
}

# Calculate LR interaction score
df.tumor.interactions$LR <- df.tumor.interactions$LIG_EXPR * 
  df.tumor.interactions$REC_EXPR
df.tumor.interactions$LR_name <- paste0(
  df.tumor.interactions$LIG, "_", 
  df.tumor.interactions$REC
)

# Rank LR interactions within each patient
df.tumor.interactions <- df.tumor.interactions %>%
  dplyr::group_by(cohort, PR, patient) %>%
  dplyr::mutate(LR_rank = rank(LR, ties.method = "random"))

# Calculate LR prevalence for GBM-source interactions
df.tumor.interactions.sum <- df.tumor.interactions %>%
  dplyr::group_by(LR_name, LIG, REC) %>%
  dplyr::summarise(
    # Combined statistics (primary + recurrent)
    LR.mean.all = mean(LR[PR %in% c("P", "R")]),
    LR.prop.p.all = mean(LR[PR %in% c("P", "R")] > 0, na.rm = TRUE),
    LR.pval.all = suppressWarnings(
      wilcox.test(LR[PR %in% c("P", "R")])$p.value
    ),
    
    # Primary GBM statistics
    LR.mean.P = mean(LR[PR %in% "P"], na.rm = TRUE),
    LR.prop.p = mean(LR[PR %in% "P"] > 0, na.rm = TRUE),
    LIG.EXPR.P = mean(LIG_EXPR[PR %in% "P"], na.rm = TRUE),
    REC.EXPR.P = mean(REC_EXPR[PR %in% "P"], na.rm = TRUE),
    
    # Recurrent GBM statistics
    LR.mean.R = mean(LR[PR %in% "R"], na.rm = TRUE),
    LR.prop.R = mean(LR[PR %in% "R"] > 0, na.rm = TRUE),
    LIG.EXPR.R = mean(LIG_EXPR[PR %in% "R"], na.rm = TRUE),
    REC.EXPR.R = mean(REC_EXPR[PR %in% "R"], na.rm = TRUE),
    
    # Change statistics
    delta = mean(LR[PR %in% "R"], na.rm = TRUE) - mean(LR[PR %in% "P"], na.rm = TRUE),
    
    # Sample counts
    n.p.active = sum(LR[PR %in% "P"] > 0, na.rm = TRUE),
    n.p.total = sum(PR %in% "P", na.rm = TRUE),
    n.r.active = sum(LR[PR %in% "R"] > 0, na.rm = TRUE),
    n.r.total = sum(PR %in% "R", na.rm = TRUE),
    
    # Statistical tests
    delta.prop.pval = suppressWarnings(
      prop.test(x = c(sum(LR[PR %in% "P"] > 0), sum(LR[PR %in% "R"] > 0)),
                n = c(sum(PR %in% "P"), sum(PR %in% "R")))[["p.value"]]
    ),
    delta.pval = suppressWarnings(
      wilcox.test(LR[PR %in% "P"], LR[PR %in% "R"])$p.value
    ),
    P.pval = suppressWarnings(wilcox.test(LR[PR %in% "P"])$p.value),
    R.pval = suppressWarnings(wilcox.test(LR[PR %in% "R"])$p.value),
    LFC = log(mean(LR[PR %in% "R"]) + 1) - log(mean(LR[PR %in% "P"]) + 1),
    LR.rank.P = mean(LR_rank[PR %in% "P"]),
    LR.rank.R = mean(LR_rank[PR %in% "R"]),
    LR.rank.pval = suppressWarnings(
      wilcox.test(LR_rank[PR %in% "P"], LR_rank[PR %in% "R"])$p.value
    ),
    .groups = "drop"
  )

# Calculate additional metrics and multiple testing correction
df.tumor.interactions.sum$prop.delta <- df.tumor.interactions.sum$LR.prop.R - 
  df.tumor.interactions.sum$LR.prop.p
df.tumor.interactions.sum$delta.fdr <- p.adjust(
  df.tumor.interactions.sum$delta.pval, method = "BH"
)
df.tumor.interactions.sum$delta.prop.fdr <- p.adjust(
  df.tumor.interactions.sum$delta.prop.pval, method = "BH"
)
df.tumor.interactions.sum$P.padj <- p.adjust(
  df.tumor.interactions.sum$P.pval, method = "BH"
)
df.tumor.interactions.sum$R.padj <- p.adjust(
  df.tumor.interactions.sum$R.pval, method = "BH"
)
df.tumor.interactions.sum$LR.padj.all <- p.adjust(
  df.tumor.interactions.sum$LR.pval.all, method = "BH"
)

message(paste0("  Calculated prevalence for ", 
               nrow(df.tumor.interactions.sum), " GBM-source LR pairs"))

# ==============================================================================
# SECTION 10: CALCULATE LIGAND PREVALENCE SCORES
# ==============================================================================

message("=== SECTION 10: Calculating ligand prevalence scores (P_L) ===")

# Create gene-to-protein mapping
gene2prot <- names(prot2gene)
names(gene2prot) <- as.character(prot2gene)

# --- OL-source ligand prevalence ---
df.olg.1 <- df.opc.sum
df.olg.1$LIG <- prot2gene[df.olg.1$cytokine]
df.olg.1$delta[df.olg.1$delta < 0] <- 0
df.olg.1 <- df.olg.1 %>% dplyr::select(c(LIG, delta))
df.olg.1$delta <- log10(df.olg.1$delta + 1)
colnames(df.olg.1) <- c("LIG", "GO1_OLG_CYTOKINE")
df.olg.1_raw <- df.olg.1

# Calculate average LR prevalence per ligand (P_L)
df.olg.2 <- df.olig.interactions.sum
df.olg.2$cytokine <- gene2prot[df.olg.2$LIG]
df.olg.2$LR.prop <- df.olg.2$LR.prop.p.all
df.olg.2 <- df.olg.2 %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(LIG, LR.prop))
colnames(df.olg.2) <- c("LIG", "GO2_OLG_LR")
df.olg.2 <- df.olg.2 %>% 
  dplyr::filter(GO2_OLG_LR != 0) %>% 
  dplyr::group_by(LIG) %>% 
  dplyr::summarise(GO2_OLG_LR = mean(GO2_OLG_LR, na.rm = TRUE))
df.olg.2_raw <- df.olg.2

# Merge OL cytokine abundance with LR prevalence
df.olg.raw <- merge(df.olg.1_raw, df.olg.2_raw, by = "LIG", all = TRUE)

# --- GBM-source ligand prevalence ---
df.go.1 <- df.g.sum
df.go.1$LIG <- prot2gene[df.go.1$cytokine]
df.go.1$delta[df.go.1$delta < 0] <- 0
df.go.1 <- df.go.1 %>% dplyr::select(c(LIG, delta))
df.go.1$delta <- log10(df.go.1$delta + 1)
colnames(df.go.1) <- c("LIG", "GO1_GBM_CYTOKINE")
df.go.1_raw <- df.go.1

# Calculate average LR prevalence per ligand
df.go.2 <- df.tumor.interactions.sum
df.go.2$cytokine <- gene2prot[df.go.2$LIG]
df.go.2$LR.prop <- df.go.2$LR.prop.p.all
df.go.2 <- df.go.2 %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(LIG, LR.prop))
colnames(df.go.2) <- c("LIG", "GO2_GBM_LR")
df.go.2 <- df.go.2 %>% 
  dplyr::filter(GO2_GBM_LR != 0) %>% 
  dplyr::group_by(LIG) %>% 
  dplyr::summarise(GO2_GBM_LR = mean(GO2_GBM_LR, na.rm = TRUE))
df.go.2_raw <- df.go.2

# Merge GBM cytokine abundance with LR prevalence
df.GO.raw <- merge(df.go.1_raw, df.go.2_raw, by = "LIG", all = TRUE)

# Harmonize gene symbols (cytokine protein names to gene symbols)
c2c <- c(
  "FGF-basic" = "FGF2", "Adiponectin" = "ADIPOQ", "VEGF" = "VEGFA", 
  "Decorin" = "DCN", "TNFa" = "TNF", "TGF-beta-1" = "TGFB1", 
  "Resistin" = "RETN", "TWEAK" = "TNFSF12", "Prolactin" = "PRL", 
  "Persephin" = "PSPN", "Noggin" = "NOG", "Cardiotrophin-1" = "CTF1",
  "Neuropoietin" = "CTF2P", "M-CSF" = "CSF1", "LIGHT" = "TNFSF14",
  "IL36RA" = "IL36RN", "IL36a" = "IL36A", "IL23" = "IL23A",
  "IL1ra" = "IL1RN", "IL1b" = "IL1B", "IL1a" = "IL1A", 
  "IL12B" = "IL12", "IL12A" = "IL12",
  "IFNk" = "IFNK", "IFNg" = "IFNG", "IFNe" = "IFNE", 
  "IFNa1" = "IFNA1", "IFNb" = "IFNB1",
  "GM-CSF" = "CSF2", "GITRL" = "TNFSF18", "G-CSF" = "CSF3", 
  "Flt3l" = "FLT3LG", "FasL" = "FASLG",
  "CD40L" = "CD40LG", "BAFF" = "TNFSF13B", "APRIL" = "TNFSF13", 
  "41BBL" = "TNFSF9"
)

# Apply name harmonization
df.olg.raw$LIG[df.olg.raw$LIG %in% names(c2c)] <- 
  c2c[df.olg.raw$LIG[df.olg.raw$LIG %in% names(c2c)]]
df.GO.raw$LIG[df.GO.raw$LIG %in% names(c2c)] <- 
  c2c[df.GO.raw$LIG[df.GO.raw$LIG %in% names(c2c)]]

# Calculate product score (cytokine abundance × LR prevalence)
df.olg.raw$prod <- df.olg.raw$GO1_OLG_CYTOKINE * df.olg.raw$GO2_OLG_LR
df.GO.raw$prod <- df.GO.raw$GO1_GBM_CYTOKINE * df.GO.raw$GO2_GBM_LR

message("  Calculated ligand prevalence scores for OL and GBM sources")

# ==============================================================================
# SECTION 11: SAVE OUTPUT FILES
# ==============================================================================

message("=== SECTION 11: Saving output files ===")

# --- Cytokine profiling results ---
write.csv(
  df.opc.sum, 
  file = file.path(config$output_dir, "cytokine_summary_OPC.csv"),
  row.names = FALSE
)

write.csv(
  df.g.sum, 
  file = file.path(config$output_dir, "cytokine_summary_GBM.csv"),
  row.names = FALSE
)

write.csv(
  df.gbm.sum, 
  file = file.path(config$output_dir, "cytokine_summary_pGBM_vs_rGBM.csv"),
  row.names = FALSE
)

# --- LR prevalence results ---
write.csv(
  df.olig.interactions.sum, 
  file = file.path(config$output_dir, "LR_prevalence_OL_source.csv"),
  row.names = FALSE
)

write.csv(
  df.tumor.interactions.sum, 
  file = file.path(config$output_dir, "LR_prevalence_GBM_source.csv"),
  row.names = FALSE
)

# --- Ligand prevalence scores ---
write.csv(
  df.olg.raw, 
  file = file.path(config$output_dir, "ligand_prevalence_OL_source.csv"),
  row.names = FALSE
)

write.csv(
  df.GO.raw, 
  file = file.path(config$output_dir, "ligand_prevalence_GBM_source.csv"),
  row.names = FALSE
)

# --- Raw interaction data (RDS format) ---
saveRDS(
  list(
    olig.interactions = df.olig.interactions,
    tumor.interactions = df.tumor.interactions,
    olig.interactions.sum = df.olig.interactions.sum,
    tumor.interactions.sum = df.tumor.interactions.sum
  ),
  file = file.path(config$output_dir, "LR_interactions_raw.rds")
)

# --- Complete results bundle ---
saveRDS(
  list(
    # Cytokine profiling
    cytokine_opc = df.opc.sum,
    cytokine_opc_top = df.opc.sum.top,
    cytokine_gbm = df.g.sum,
    cytokine_gbm_top = df.g.sum.top,
    cytokine_pr = df.gbm.sum,
    
    # LR prevalence
    lr_olig_source = df.olig.interactions.sum,
    lr_tumor_source = df.tumor.interactions.sum,
    
    # Ligand prevalence
    ligand_olig = df.olg.raw,
    ligand_gbm = df.GO.raw,
    
    # Mappings
    prot2gene = prot2gene,
    gene2prot = gene2prot,
    
    # LR database subset
    lr_database = df.LR.subset
  ),
  file = file.path(config$output_dir, "PR14_complete_results.rds")
)

message(paste0("  Saved results to: ", config$output_dir))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
