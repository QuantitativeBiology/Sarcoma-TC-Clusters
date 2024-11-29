# Load Libraries
library(ConsensusClusterPlus)
library(readxl)
library(ComplexHeatmap)
library(dplyr)
library(tidyverse)
library(reshape2)
library(GGally)
library(circlize)

# Define output file path
output_file <- "RESULTS/Consensus_Clusters.csv"

# Check if the output file exists
if (!file.exists(output_file)) {
  # Read normalized counts
  voom <- read.csv("FILES/normalized_counts.csv", row.names = 1)
  colnames(voom) <- sub("^X", "", colnames(voom))
  
  # Define title for the consensus results
  title <- "RESULTS/ConsensusResults"
  
  # Calculate MAD for each row
  row_sds <- apply(voom, MARGIN = 1, mad)
  
  # Sort by descending MAD
  row_sds <- row_sds[order(row_sds, decreasing = TRUE)]
  
  # Select top percentage of rows
  top_percentage <- 0.55
  barplot(row_sds, xlab = "Genes", ylab = "Standard Deviation", las = 2)
  num_rows_to_keep <- ceiling(top_percentage * length(row_sds))
  row_sds <- row_sds[1:num_rows_to_keep]
  
  # Filter voom data to keep selected rows
  voom <- voom[names(row_sds), ]
  
  # Centre rows by subtracting their median
  voom <- sweep(voom, 1, apply(voom, 1, median, na.rm = TRUE))
  
  # Perform consensus clustering
  results <- ConsensusClusterPlus(
    data.matrix(voom),
    maxK = 10,
    reps = 50,
    pItem = 0.8,
    pFeature = 1,
    title = title,
    clusterAlg = "hc",
    distance = "pearson",
    seed = 1262118388.71279,
    plot = "png"
  )
  
  # Calculate item-consensus matrix
  icl <- calcICL(results, title = title, plot = "png")
  
  # Extract the consensus classes
  final <- results[[4]][["consensusClass"]]
  values <- data.frame(final)
  
  # Save results to CSV
  write.csv(values, output_file)
  
  print("Consensus clusters calculated and saved.")
} else {
  print("Output file already exists. No action taken.")
}