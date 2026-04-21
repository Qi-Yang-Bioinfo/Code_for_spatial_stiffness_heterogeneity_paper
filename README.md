# Confinement creates tissue mechanics zonation to spatially pattern and accelerate solid tumor growth

**Repository:** [Code_for_spatial_stiffness_heterogeneity_paper](https://github.com/Qi-Yang-Bioinfo/Code_for_spatial_stiffness_heterogeneity_paper)

This repository contains the standalone R scripts and source code used for the single-cell/nuclear RNA-sequencing (scRNA-seq/snRNA-seq) and bulk transcriptome data analysis in the manuscript *"Confinement creates tissue mechanics zonation to spatially pattern and accelerate solid tumor growth"*. The pipeline covers raw count pre-processing, doublet removal, sample integration, clustering, cell-type annotation, pathway scoring for both in-house medulloblastoma models and public human datasets.

---

## 1. System Requirements

### Hardware Requirements
The basic scripts can run on a standard desktop computer. However, the full scRNA-seq integration, clustering, CytoTRACE analyses require substantial computational resources.
* **Optimal system:** A High-Performance Computing (HPC) node or workstation with a minimum of 64 GB RAM and 16+ CPU cores. 
* **Non-standard hardware:** None required.

### Software Requirements & OS
The scripts have been written in **R** and tested on **Linux (Ubuntu 22.04.5 LTS)**, but they are compatible with macOS and Windows environments supporting R.

**Tested Versions:**
* **R version:** 4.4.0 or higher
* **Python version:** 3.8 or higher (Required for CytoTRACE via `reticulate`)

**Major R Package Dependencies:**
* `Seurat` (v5.0+)
* `DoubletFinder` (v2.0.4)
* `CytoTRACE` (v0.3.3, requires Python environment setup)

---

## 2. Installation Guide

### Instructions

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Qi-Yang-Bioinfo/Code_for_spatial_stiffness_heterogeneity_paper.git
   cd Code_for_spatial_stiffness_heterogeneity_paper

## 3. Datasets used in this manuscript

* **In-house snRNA-seq dataset:** Raw data and processed data can be accessed under [GSE315300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315300).
* **Human MB bulk transcriptome dataset:** [GSE124814](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124814)
* **Human MB scRNA-seq transcriptome dataset:** [GSE156053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156053)
* **Human basal cell carcinoma scRNA-seq dataset:** [GSE141526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141526)

Code and examples for processing these datasets can be found in the corresponding scripts.

