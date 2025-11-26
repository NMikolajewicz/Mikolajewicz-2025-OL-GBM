# Mikolajewicz-2025-OL-GBM

**Reactive oligodendrocytes promote glioblastoma progression through CCL5/CCR5-mediated glioma stem cell maintenance**

[![DOI](https://img.shields.io/badge/DOI-10.xxxx%2Fxxxxx-blue)](https://doi.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the reproducible R code for single-cell RNA sequencing analyses presented in Mikolajewicz et al. 2025.

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Data Acquisition](#data-acquisition)
- [Usage](#usage)
- [Analysis Workflow](#analysis-workflow)
- [Output Files](#output-files)
- [Configuration](#configuration)
- [Runtime Estimates](#runtime-estimates)
- [Citation](#citation)
- [License](#license)

## Overview

This study investigates the role of oligodendrocytes (OLs) in glioblastoma (GBM) progression using:

- **Single-cell RNA sequencing** of primary (pGBM) and recurrent (rGBM) tumors
- **Cytokine profiling** (96-plex immunoassays)
- **Ligand-receptor interaction analysis**
- **Spatial transcriptomics** validation

Key findings:
- OLs are enriched in the GBM tumor microenvironment, particularly at recurrence
- GBM cells recruit OLs via CX3CL1/CX3CR1 signaling
- Interferon-driven reactive OL states secrete pro-tumorigenic cytokines (notably CCL5)
- CCL5/CCR5 signaling promotes glioma stem cell maintenance

## Repository Structure

```
Mikolajewicz-2025-OL-GBM/
├── R/                              # Core analysis functions
│   ├── utils.R                     # Utility functions (I/O, logging)
│   ├── qc.R                        # Quality control and preprocessing
│   ├── integration.R               # Multi-sample integration (rPCA, BBKNN)
│   ├── cnv.R                       # CNV-based tumor calling
│   └── differential_expression.R   # DEG and abundance analysis
│
├── scripts/                        # Analysis run scripts
│   └── 01_GBM_annotation.R         # Main GBM analysis pipeline
│
├── config/
│   └── config.yml                  # Centralized parameters and paths
│
├── data/
│   ├── raw/                        # Input Seurat objects (not tracked)
│   ├── processed/                  # Processed outputs
│   └── reference/                  # Gene sets, LR databases
│
├── results/
│   ├── figures/                    # Generated figures
│   └── tables/                     # Generated tables
│
├── renv.lock                       # Package versions for reproducibility
└── README.md                       # This file
```

## Installation

### Prerequisites

- R ≥ 4.2.0
- Python ≥ 3.8 (for BBKNN integration, optional)

### R Package Installation

```r
# Install renv for reproducible environment
install.packages("renv")

# Clone repository and restore environment
# From terminal:
# git clone https://github.com/NMikolajewicz/Mikolajewicz-2025-OL-GBM.git
# cd Mikolajewicz-2025-OL-GBM

# Restore packages from renv.lock
renv::restore()
```

### Required R Packages

The following packages are used (versions in `renv.lock`):

**Core Analysis:**
- Seurat (≥5.0.0)
- sctransform (≥0.4.0)
- glmGamPoi
- presto
- mclust
- harmony

**Data Manipulation:**
- dplyr
- tidyr
- Matrix

**Visualization:**
- ggplot2
- cowplot
- pheatmap

**Integration:**
- reticulate (for BBKNN)
- msigdbr

### Python Dependencies (Optional)

For BBKNN batch correction:

```bash
pip install bbknn scanpy anndata
```

## Data Acquisition

### In-house GBM Cohort

Sci-RNA-seq3 data generated in this study is available on Figshare:
- DOI: [10.6084/m9.figshare.25917628.v1](https://doi.org/10.6084/m9.figshare.25917628.v1)

Download and place files in `data/raw/`:

```bash
# Download from Figshare
wget -P data/raw/ <figshare_url>/p9_hs_pr_gbm.rds
wget -P data/raw/ <figshare_url>/p14_hs_pr_gbm.rds
wget -P data/raw/ <figshare_url>/p16_hs_pr_gbm.rds
```

### External Cohorts

The following public datasets are used for validation:

| Cohort | Source | Accession |
|--------|--------|-----------|
| Abdelfattah | GEO | GSE182109 |
| Wang | GEO | GSE131928 |
| Ravi (spatial) | Datadryad | 10.5061/dryad.h70rxwdmj |

## Usage

### Quick Start

```r
# Set working directory to repository root
setwd("path/to/Mikolajewicz-2025-OL-GBM")

# Run main analysis pipeline
source("scripts/01_GBM_annotation.R")
```

### Step-by-Step Execution

```r
# 1. Load configuration and functions
source("R/utils.R")
config <- load_config("config/config.yml")
source("R/qc.R")
source("R/integration.R")
source("R/cnv.R")
source("R/differential_expression.R")

# 2. Load data
so <- readRDS("data/raw/p9_hs_pr_gbm.rds")

# 3. Preprocess
so <- preprocess_sample(so, config, sample_name = "p9")

# 4. Continue with analysis...
```

## Analysis Workflow

The analysis reproduces the following manuscript figures:

### Figure 1: OL Enrichment in GBM
- **1A-B**: Cell type fractional abundance (pGBM vs rGBM)
- **1C-E**: Spatial distribution of OLs

### Figure 2: CX3CL1/CX3CR1 Signaling
- **2A-B**: Cytokine profiling and LR interactions

### Figure 3: OL Meta-Atlas
- **3A-D**: Pan-disease OL metaprograms (NMF analysis)

### Key Analysis Steps

1. **Quality Control**
   - Gene/cell filtering: 200-9000 genes, <60% mito
   - SCTransform normalization (v2, glmGamPoi)

2. **Integration**
   - Reciprocal PCA (rPCA) across samples
   - BBKNN for batch visualization

3. **Tumor Classification**
   - PCA-based CNV analysis (chromosome arm scoring)
   - Microglia-based reference correction
   - GMM-based classification with posterior probability

4. **Cell Type Annotation**
   - Module scoring for canonical markers
   - Over-clustering validation

5. **Differential Expression**
   - Wilcoxon tests via presto
   - Abundance comparisons (paired, unpaired)

## Output Files

### Tables

| File | Description |
|------|-------------|
| `TableS_tumor_classification_summary.csv` | Tumor/normal cell counts per sample |
| `TableS_celltype_summary.csv` | Cell type counts and percentages |
| `TableS_abundance_PR_comparison.csv` | Primary vs recurrent abundance stats |
| `TableS_normal_celltype_markers.csv` | DE markers for non-tumor cell types |
| `TableS_tumor_cluster_markers.csv` | DE markers for tumor clusters |

### Figures

| File | Description |
|------|-------------|
| `Fig1A_fractional_abundance.png` | Cell type abundances (pGBM vs rGBM) |
| `FigS1_UMAP_overview.png` | UMAP colored by cell type, sample, tumor status |
| `FigS1_marker_dotplot.png` | Marker gene expression dot plot |

### Processed Data

| File | Description |
|------|-------------|
| `so_integrated.rds` | Integrated Seurat object with all annotations |
| `so_tumor.rds` | Tumor cell subset |
| `so_normal.rds` | Non-tumor cell subset |

## Configuration

All parameters are centralized in `config/config.yml`:

```yaml
# Key parameters (matching manuscript methods)
qc:
  min_genes: 200
  max_genes: 9000
  max_mito_percent: 60

dim_reduction:
  umap:
    n_neighbors: 50
    min_dist: 0.1
    metric: "cosine"

cnv:
  posterior_prob_threshold: 0.5
  reference_clusters: [3, 5, 16, 13, 15]  # Microglia
```

Modify `config.yml` to adjust parameters or paths.

## Runtime Estimates

| Analysis Step | Cells | Time (approx.) |
|--------------|-------|----------------|
| Preprocessing (per sample) | ~10,000 | 5-10 min |
| Integration (all samples) | ~50,000 | 30-60 min |
| CNV analysis | ~50,000 | 15-30 min |
| Differential expression | ~50,000 | 5-10 min |
| **Total** | **~50,000** | **~2 hours** |

Tested on: 16-core CPU, 64GB RAM

## Reproducibility

To ensure exact numerical reproducibility:

1. Use the provided `renv.lock` to restore package versions
2. Random seed is set via `config$reproducibility$seed` (default: 123)
3. All key parameters are documented in `config.yml`

```r
# Verify reproducibility
set.seed(123)
sessionInfo()
```

## Citation

If you use this code or data, please cite:

```
Mikolajewicz N, Zhai K, et al. (2025). Reactive oligodendrocytes promote 
glioblastoma progression through CCL5/CCR5-mediated glioma stem cell 
maintenance. [Journal]. DOI: [pending]
```

## Contact

- **Lead Contact**: Jason Moffat (jason.moffat@sickkids.ca)
- **Code Issues**: [GitHub Issues](https://github.com/NMikolajewicz/Mikolajewicz-2025-OL-GBM/issues)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Troubleshooting

### Common Issues

**1. BBKNN not available**
```
Warning: BBKNN not available, falling back to Harmony...
```
Install Python dependencies: `pip install bbknn scanpy`

**2. Memory errors during integration**
```r
# Reduce memory usage:
config$normalization$sctransform$conserve_memory <- TRUE
options(future.globals.maxSize = 16000 * 1024^2)  # 16GB
```

**3. Missing genes in module scoring**
```
Warning: Only 60% of requested genes found
```
This is expected for some marker sets. The analysis continues with available genes.

### Logs

All operations are logged with timestamps. Check console output or redirect to file:

```r
sink("analysis_log.txt", split = TRUE)
source("scripts/01_GBM_annotation.R")
sink()
```
