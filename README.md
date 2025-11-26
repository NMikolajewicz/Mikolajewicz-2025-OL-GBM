# README

# Reactive Oligodendrocytes in Glioblastoma Progression

Mikolajewicz et al. 2025, _Neuron_

## Overview

This repository contains analysis code for studying oligodendrocyte-glioblastoma interactions through CCL5/CCR5-mediated signaling. The work integrates single-cell RNA sequencing, cytokine profiling, and ligand-receptor interaction analysis across primary and recurrent GBM samples.

## Repository Structure

```
├── GBM_scRNAseq_analysis/     # scRNA-seq preprocessing, QC, integration, and tumor calling
├── OL_meta_atlas/             # Oligodendrocyte meta-atlas with NMF metaprogram discovery
├── Cytokines/                 # Cytokine profiling and ligand-receptor prevalence analysis
```

## Datasets

| Dataset | DOI |
|---------|-----|
| snRNA-seq Primary-Recurrent GBM (Mikolajewicz Cohort) | [10.6084/m9.figshare.25917628](https://figshare.com/articles/dataset/snRNA-seq_Primary-Recurrent_GBM_Mikolajewicz_Cohort_/25917628) |
| Oligodendrocyte Meta-Atlas | [10.6084/m9.figshare.30702419](https://figshare.com/articles/dataset/scRNA-seq_Oligodendrocyte-Lineage_Meta-Atlas/30702419) |

## Requirements

R with Seurat, sctransform, dplyr, tidyr, ggplot2, and NMF packages. BBKNN batch correction requires Python ≥ 3.7 with scanpy and bbknn. See individual module READMEs for complete dependencies.

## Citation

Mikolajewicz, N., Zhai, K., Puri, A., Miletic, P., Tatari, N., Wei, J., Savage, N., Huang, Z., Huang, Q., Lee, S. Y., Ahmadnejad, M. J., Nguyen, R., Chen, D., Korman, T., Mobilio, D., Topley, M., Lu, J. Q., Voisin, M. R., Zador, Z., Chafe, S. C., Venugopal, C., Brown, K. R., Zadeh, G., Han, H., Muffat, J., Bao, S., Singh, S. K., & Moffat, J. (2025). Reactive oligodendrocytes promote glioblastoma progression through CCL5/CCR5-mediated glioma stem cell maintenance. _Neuron_.

## Contact

Nicholas Mikolajewicz — n.mikolajewicz@utoronto.ca
