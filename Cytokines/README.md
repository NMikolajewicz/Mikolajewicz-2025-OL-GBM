# OL-GBM Cytokine Profiling and Ligand-Receptor Analysis

## Overview

This repository contains the consolidated R script for cytokine profiling and ligand-receptor (LR) interaction prevalence analysis performed in the OL-GBM study:

1. **Quantify cytokine abundance** in conditioned media from OL and GBM cells
2. **Compare cytokine profiles** between primary and recurrent GBM samples
3. **Generate ligand prevalence scores** combining cytokine abundance with LR interaction frequency

The main script (`ligand_receptor_analysis.R`) performs all data loading, preprocessing, statistical analysis, and output generation in a single reproducible workflow.

## Repository Structure

```
├── ligand_receptor_analysis.R           # Main R script
├── README.md                            # This file
```

## Requirements

### Required R Packages

**Core packages** (available from CRAN):
```r
# Data manipulation and analysis
Seurat >= 4.0.0
scMiko >= 1.0.0
sctransform >= 0.3.5
dplyr >= 1.0.0
tidyr >= 1.0.0

# Data import
readxl >= 1.4.0

# Visualization
ggplot2 >= 3.3.0
ggrepel >= 0.9.0
RColorBrewer >= 1.1-2
circlize >= 0.4.16
reshape2 >= 1.4.4
cowplot >= 1.1.0
```

### Installation Instructions
```r
# Install CRAN packages
install.packages(c(
  "Seurat", "sctransform", "dplyr", "tidyr", "readxl",
  "ggplot2", "ggrepel", "RColorBrewer", "reshape2", 
  "circlize", "cowplot"
))

# Install scMiko package from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github(repo = "NMikolajewicz/scMiko")
```

### Key scMiko Functions Used
The scMiko package provides essential helper functions:
- `LR.db` - Curated ligand-receptor interaction database
- `theme_miko()`, `scale_color_miko()`, `scale_fill_miko()` - Custom plotting themes
- `p2z()`, `z2p()` - Statistical conversion utilities
- `categoricalColPal()` - Color palette generation

## Data Inputs

### 1. Cytokine Profiling Data
**File**: `data/cytokine_analysis_280224.xlsx`  
**Format**: Excel workbook with "clean_data" sheet  
**Contents**:
- 96-plex immunoassay results from conditioned media experiments
- Columns: Sample metadata (ID, Sample1, Sample2, Sample3) + 96 cytokine measurements
- Samples include:
  - OPC conditioned media (OPC CM)
  - GBM conditioned media (GBM CM)
  - Primary vs. recurrent GBM comparisons

**Notes**:
- "OOR <" values (out of range low) are automatically converted to 0
- Measurements are in pg/mL

### 2. Single-Cell RNA-seq Data
**File**: `data/PR4_preprocessed_tumor_and_other_miko_wang_abdel_140623.rds`  
**Format**: R data file containing preprocessed Seurat objects  
**Contents**:
- Integrated scRNA-seq data from multiple GBM cohorts:
  - Mikolajewicz cohort (primary source)
  - Abdelfattah cohort (validation)
  - Wang cohort (validation)
- Cell type annotations (tumor cells, oligodendrocytes, immune cells, etc.)
- Patient metadata including primary/recurrent status
- Normalized expression matrices

### Configuration

Modify file paths in the script's configuration section:

```r
config <- list(
  # Update these paths according to your directory structure
  cytokine_data_path = "data/cytokine_analysis_280224.xlsx",
  seurat_data_path = "data/PR4_preprocessed_tumor_and_other_miko_wang_abdel_140623.rds",
  output_dir = "outputs/PR14_results/"
)
```
