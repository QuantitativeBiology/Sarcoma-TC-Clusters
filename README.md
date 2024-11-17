
# Sarcoma-TC-Clusters

This repository contains scripts and data for the analysis for the publication "Machine learning-based analysis of genomic and transcriptomic data unveils sarcoma clusters with superlative prognostic and predictive value". 

---

## Introduction

Sarcomas are a heterogeneous group of cancers with distinct genetic profiles. This repository explores the identification and validation of transcriptomic clusters in sarcoma. 

Note - Proprietary dataset from IPO used can be provided upon request

The repository leverages TCGA-SARC data and integrates external validation datasets to ensure reproducibility and robustness of the results.

---

## Folder Structure

```plaintext
├── FILES
│   ├── Proprietary files can be sent upon request
├── RESULTS
│   ├── TCGA-TC-Classified.csv # TCGA-SARC Patients Classified with TC
│   ├── Consensus_Clusters_Genes.txt # Genes discovered from Proprietary Cohort that define clusters
│   ├── clinical_cinsarc_paper.csv # TCGA-SARC classified with CINSARC
│   ├── Sarculator_TCGA.csv # TCGA-SARC Classified with SARCULATOR
│   ├── TCGA-SARC.csv # TCGA-SARC Gene Expression Matrix
│   └── Sarculator_Results.xlsx # Proprietary patients classified with TC
├── *.R                               # Main R scripts for analysis
├── .gitignore                        # Ignore configuration
└── README.md                         # Project documentation
```

---

## Dependencies

To reproduce the analyses in this repository, the following R packages are required:

- `tidyverse`
- `limma`
- `edgeR`
- `ComplexHeatmap`
- `EnhancedVolcano`
- `clusterProfiler`
- `msigdbr`
- `fgsea`
- `survival`
- `survminer`

Install them in R using the following command:
```R
install.packages(c("tidyverse", "limma", "edgeR", "ComplexHeatmap", "EnhancedVolcano", 
                   "clusterProfiler", "msigdbr", "fgsea", "survival", "survminer"))
```

---

## Usage

Scripts must be run in ascendant order

### 0. Data Preprocessing
The script `0.clinical_data_preprocessing.R` processes and cleans the clinical data and RNA expression data.

### 1. Consensus Clustering
Run `1.consensus_cluster.R` to identify transcriptomic clusters and extract cluster-specific genes.

### 2. Survival Analysis
Use `2.survival_analysis.R` to evaluate the clinical relevance of the identified clusters using Kaplan-Meier survival plots and Cox proportional hazard models.

### 3. TCGA-SARC Validation
Scripts in the `3.*.R` series validate the clusters using TCGA-SARC, integrating and comparing with nomogram tools - Sarculator and gene expression signatures - CINSARC.


## Citation

To be filled

---

## License

This repository is licensed under the MIT License. See `LICENSE` for details.
