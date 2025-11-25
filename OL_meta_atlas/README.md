# Oligodendrocyte Meta-Atlas Analysis

## Overview

This repository contains the R code for generating an oligodendrocyte meta-atlas by integrating single-cell/nucleus RNA sequencing data from multiple studies spanning healthy brain, neurological diseases, and brain tumors. The analysis identifies robust transcriptional programs (metaprograms) using Non-negative Matrix Factorization (NMF).

## Repository Structure

```
├── OL_meta_atlas_analysis.R             # Main R script (preprocessing + integration + NMF)
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
| Study | GEO/Source | Journal |
|-------|-----------|----------|
| Franjic et al., 2022 | GSE186538 | Neuron |
| Habib et al., 2017 | - | Nature Methods |
| Kanton et al., 2019 | E-MTAB-8230 | Nature |
| Khrameeva et al., 2020 | GSE127774 | Genome Research |
| Hodge et al., 2019 | Allen Brain Atlas | Nature |
| Bakken et al., 2021 | SCR_016152 | Nature |

### Multiple Sclerosis
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Wheeler et al., 2020 | GSE130119 | Nature |
| Jäkel et al., 2019 | GSE118257 | Nature |
| Schirmer et al., 2019 | Cell Browser (ms) | Nature |

### Brain Tumors
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Sun et al., 2022 | GSE202371 | Clinical and Translational Medicine |
| Kim et al., 2020 | GSE131907 | Nature Communications |
| Biermann et al., 2022 | GSE185386 | Cell |
| Heming et al., 2022 | GSE203552 | Genome Medicine |

### Neurodegenerative Diseases
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Lau et al., 2020 | GSE157827 | PNAS |
| Smajić et al., 2022 | GSE157783 | Brain |

### Developmental
| Study | GEO/Source | Citation |
|-------|-----------|----------|
| Yu et al., 2021 | GSE165388 | Nature Neuroscience |
| van Bruggen et al., 2022 | Zenodo | Developmental Cell |
| Cao et al., 2020 | - | Science |
| Bhaduri et al., 2021 | SCR_002001 | Nature |
| Aldinger et al., 2021 | Cell Browser (cbl-dev) | Nature Neuroscience |

### Glioblastoma (Primary Analysis)
| Study | Source | Citation |
|-------|--------|----------|
| Mikolajewicz et al., 2025 | This study | - |
| Abdelfattah et al., 2022 | GSE182109 | Nature Communications |
| Wang et al., 2021 | GSE131928 | Genome Biology |

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

## Contact

For questions about the analysis, please contact [n.mikolajewicz@utoronto.ca].

