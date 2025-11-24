# Oligodendrocyte Meta-Atlas Analysis

## Overview

This repository contains the R code for generating an oligodendrocyte meta-atlas by integrating single-cell/nucleus RNA sequencing data from multiple studies spanning healthy brain, neurological diseases, and brain tumors. The analysis identifies robust transcriptional programs (metaprograms) using Non-negative Matrix Factorization (NMF).

## Associated Publication

Mikolajewicz N, et al. *[Manuscript Title]*

## Repository Structure

```
├── OL_meta_atlas_analysis.Rmd          # Main R Notebook (integration + NMF)
├── OL_preprocessing_individual_studies.R # Preprocessing for each study
├── README.md                            # This file
```

## Requirements

### R Version
- R >= 4.0.0

### Python (for BBKNN batch correction)
- Python >= 3.7
- bbknn (`pip install bbknn`)
- scanpy (`pip install scanpy`)

### Required R Packages

```r
# Core scRNA-seq analysis
install.packages(c("Seurat", "sctransform", "glmGamPoi"))

# NMF analysis
install.packages(c("NNLM", "NMF"))

# Data manipulation
install.packages(c("dplyr", "tidyr", "reshape2", "stringr"))

# Visualization  
install.packages(c("ggplot2", "RColorBrewer", "viridis", "cowplot", "pheatmap", "scattermore"))

# Parallel processing
install.packages(c("foreach", "parallel", "doParallel", "pbapply"))

# Additional
install.packages(c("Matrix", "DiffCorr"))

# Python interface
install.packages("reticulate")

# Optional: Alternative batch correction (if BBKNN unavailable)
install.packages("harmony")

# For gene annotation (Khrameeva study)
BiocManager::install("biomaRt")
```

### Python Setup for BBKNN

```bash
# Create conda environment (recommended)
conda create -n bbknn python=3.9
conda activate bbknn
pip install bbknn scanpy numpy

# Configure reticulate in R
reticulate::use_condaenv("bbknn")
```

## Data Availability

Input data consists of preprocessed Seurat objects containing oligodendrocyte cells extracted from the following studies:

### Healthy Brain
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Franjic et al., 2022 | GSE178265 | Science 2022 |
| Habib et al., 2017 | GSE97930 | Nature Methods 2017 |
| Kanton et al., 2019 | - | Nature 2019 |
| Khrameeva et al., 2020 | GSE127774 | Genome Research 2020 |
| Hodge et al., 2019 | Allen Brain Atlas | Nature 2019 |
| Bakken et al., 2021 | - | Nature 2021 |

### Multiple Sclerosis
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Wheeler et al., 2020 | GSE138852 | Nature 2020 |
| Jäkel et al., 2019 | GSE118257 | Nature 2019 |
| Schirmer et al., 2019 | GSE118257 | Nature 2019 |

### Brain Tumors
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Sun et al., 2022 | - | Brain metastasis/Glioma |
| Kim et al., 2020 | GSE131907 | Brain metastasis |
| Biermann et al., 2022 | - | Brain metastasis |
| Heming et al., 2022 | GSE203187 | PCNSL |

### Neurodegenerative Diseases
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Lau et al., 2020 | GSE138852 | Alzheimer's Disease |
| Smajić et al., 2022 | - | Parkinson's Disease |

### Developmental
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Yu et al., 2021 | - | Dev brain |
| van Bruggen et al., 2022 | - | Dev brain |
| Cao et al., 2020 | - | Human cell atlas |
| Bhaduri et al., 2021 | - | Dev brain |
| Aldinger et al., 2021 | - | Cerebellum |

### Glioblastoma (Primary Analysis)
| Study | Source | Citation |
|-------|--------|----------|
| Mikolajewicz et al., 2024 | This study | - |
| Abdelfattah et al., 2022 | - | GBM |
| Wang et al., 2022 | - | GBM |

## Usage

### Step 1: Preprocessing Individual Studies

First, run the preprocessing script to extract oligodendrocytes from each study:

```r
# Set your data directory
base_dir <- "path/to/raw_data/"

# Source the preprocessing functions
source("OL_preprocessing_individual_studies.R")

# Process all studies
all_samples <- process_all_studies(base_dir, reprocess = TRUE)
```

### Step 2: Run Main Analysis

Open the R Notebook and run all chunks:

```r
# In RStudio
rmarkdown::render("OL_meta_atlas_analysis.Rmd")
```

Or run section by section interactively.

### Step 3: Configure Paths

Edit the `config` list in the notebook:

```r
config <- list(
  data_dir = "path/to/preprocessed_data/",
  output_dir = "path/to/output/",
  # ... other parameters
)
```

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────┐
│                  PREPROCESSING PHASE                        │
├─────────────────────────────────────────────────────────────┤
│  1. Load raw data from each study (23 studies)             │
│  2. QC filtering (mito < 10%, 200 < genes < 9000)         │
│  3. Identify oligodendrocyte populations                   │
│  4. SCTransform normalization per sample                   │
│  5. Save preprocessed Seurat objects                       │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                   INTEGRATION PHASE                         │
├─────────────────────────────────────────────────────────────┤
│  6. Merge all samples into single object                   │
│  7. Standard normalization on merged data                  │
│  8. PCA dimensionality reduction                           │
│  9. BBKNN batch correction (Python via reticulate)         │
│  10. UMAP visualization                                    │
│  11. SCTransform for NMF input                             │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                     NMF PHASE                               │
├─────────────────────────────────────────────────────────────┤
│  12. Per-sample NMF (k = 2-15, parallel processing)        │
│  13. Extract top 50 genes per NMF component                │
│  14. Within-sample robustness (Jaccard > 0.7)              │
│  15. Cross-sample robustness (≥3 samples)                  │
│  16. Hierarchical clustering of robust programs            │
│  17. Gene prevalence filtering (> 0.3)                     │
│  18. Final 8 metaprograms defined                          │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                    OUTPUT PHASE                             │
├─────────────────────────────────────────────────────────────┤
│  19. Module scoring for all cells                          │
│  20. Export Table S5 (gene programs)                       │
│  21. Generate visualizations                               │
│  22. Save final Seurat object                              │
└─────────────────────────────────────────────────────────────┘
```

## Output Files

| File | Description |
|------|-------------|
| `TableS5_OL_metaprograms.csv` | Gene lists for each oligodendrocyte metaprogram |
| `seurat_oligo_integrated.rds` | Integrated Seurat object (before module scoring) |
| `seurat_oligo_atlas_final.rds` | Final Seurat object with module scores |
| `oligo_nmf_programs.rds` | R list object of gene programs |
| `nmf_results_per_sample.rds` | Per-sample NMF results (for reproducibility) |
| `nmf_results_final.rds` | Final NMF results including Jaccard matrices |
| `cell2sample_mapping.rds` | Cell barcode to sample name mapping |
| `umap_by_type.pdf` | UMAP colored by condition |
| `umap_by_study.pdf` | UMAP colored by study |
| `umap_module_scores.pdf` | UMAP colored by all module scores |
| `program_jaccard_heatmap.pdf` | Program similarity heatmap |
| `program_correlations.pdf` | Program correlation heatmap |

## Identified Gene Programs

The analysis identifies 8 oligodendrocyte transcriptional programs:

| Program | Annotation | Key Markers | Description |
|---------|------------|-------------|-------------|
| O1-Neuro-I | Neuronal-like | ACTN2, SLC5A11 | Neuronal-like signature |
| O2-Reactive-I | Immune-reactive | CD74, S100A11 | Complement, phagocytosis |
| O3-Cycling | Proliferation | TOP2A, NUF2 | Cell cycle, proliferation |
| O4-OPC | OPC | PTPRZ1, SOX6 | Oligodendrocyte progenitors (PDGFRA+) |
| O5-Reactive-II | Interferon | IFI6, B2M | Interferon-responsive |
| O6-Myelin | Myelination | PLP1, MBP | Mature myelinating OLs |
| O7-Stress | Stress response | FOS, UBC | Heat shock, immediate early genes |
| O8-Neuro-II | Neuronal-like II | SNAP25, NRG3 | Synapse-associated |

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_cells_per_sample` | 50 | Minimum cells to retain a sample |
| `mito_threshold` | 10% | Maximum mitochondrial content |
| `min_features` | 200 | Minimum genes per cell |
| `max_features` | 9000 | Maximum genes per cell |
| `n_variable_features` | 2000 | Variable features for SCT |
| `nmf_k_range` | 2-15 | Range of k values for NMF |
| `top_n_genes` | 50 | Top genes per NMF component |
| `intra_threshold` | 0.7 | Within-sample Jaccard threshold |
| `inter_threshold` | 0.2 | Between-sample Jaccard threshold |
| `min_samples_robust` | 3 | Min samples for cross-sample robustness |
| `prevalence_threshold` | 0.3 | Gene prevalence for final programs |
| `n_final_programs` | 8 | Number of final gene programs |

## Troubleshooting

### BBKNN Installation Issues

If BBKNN is not available, the script automatically falls back to Harmony:

```r
# Check if BBKNN is available
reticulate::py_module_available("bbknn")

# Install BBKNN in Python
system("pip install bbknn scanpy")
```

### Memory Issues

For large datasets, use conservative memory settings:

```r
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
gc()
```

### Parallelization

Adjust number of cores based on your system:

```r
config$n_cores <- parallel::detectCores() - 2
```

## Citation

If you use this code or data, please cite:

```
[Citation information to be added upon publication]
```

## Contact

For questions about the analysis, please contact [corresponding author email].

## License

[License information]
