
# Function to install a CRAN package if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# List of required CRAN packages
cran_packages <- c(
  "readxl", "dplyr", "tidyverse", "reshape2", "GGally",
  "circlize", "forestploter", "gplots", "PerformanceAnalytics"
)

# Install any missing CRAN packages
invisible(lapply(cran_packages, install_if_missing))

# Ensure BiocManager is installed to manage Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required Bioconductor packages
bioc_packages <- c(
  "limma", "edgeR", "ComplexHeatmap", "EnhancedVolcano", 
  "clusterProfiler", "fgsea", "msigdbr", "corto", 
  "survival", "survminer", "ConsensusClusterPlus"
)

# Install any missing Bioconductor packages
invisible(lapply(bioc_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
}))
