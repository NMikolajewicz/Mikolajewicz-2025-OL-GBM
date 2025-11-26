# Mikolajewicz-2025-OL-GBM

[![DOI](https://img.shields.io/badge/DOI-10.6084/m9.figshare.25917628-blue)](https://figshare.com/articles/dataset/snRNA-seq_Primary-Recurrent_GBM_Mikolajewicz_Cohort_/25917628)

## Overview

R scripts used to run preprocessing and annotation analysis of GBM scRNA-seq datasets. 

## Repository Structure

```
Mikolajewicz-2025-OL-GBM/
├── GBM_scRNAseq_analysis/          # Core analysis functions
│   ├── utils.R                     # Utility functions (I/O, logging)
│   ├── qc.R                        # Quality control and preprocessing
│   ├── integration.R               # Multi-sample integration (rPCA, BBKNN)
│   ├── cnv.R                       # CNV-based tumor calling
│   └── differential_expression.R   # DEG and abundance analysis
│   └── plotting.R                  # Some plotting utility functions
│   └── main.R                      # Main GBM analysis pipeline
│   └── config.yml                  # Centralized parameters and paths
│   └── README.md                   # this file
├── ...                             # other analysis directories
```

## Analysis Workflow

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

## Configuration

All parameters are centralized in `config/config.yml`:
Modify `config.yml` to adjust parameters or paths.

## Contact

- **Authoer**: Nicholas Mikolajewicz (n.mikolajewicz@utoronto.ca)
