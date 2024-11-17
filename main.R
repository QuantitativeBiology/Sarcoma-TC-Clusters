# Main script to execute all steps sequentially

# Define the list of scripts in execution order
scripts <- c(
  "0.clinical_data_preprocessing.R",
  "0.rna_normalization.R",
  "1.0.consensus_cluster.R",
  "1.1.consensus_cluster_genes_extraction.R",
  "1.2.cluster_consensus_plots.R",
  "1.3.alterations_per_cluster.R",
  "2.survival_analysis.R",
  "3.0.TCGA-Validation.R",
  "3.1.TCGA-Validation-Sarculator-Cinsarc.R"
)

# Sequentially source each script
for (script in scripts) {
  cat("Running script:", script, "\n")
  source(script)
  cat("Finished running script:", script, "\n\n")
}

cat("All scripts executed successfully.\n")