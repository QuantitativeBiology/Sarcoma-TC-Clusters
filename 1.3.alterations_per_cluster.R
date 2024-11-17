library(tidyverse)
library(reshape2)
library(gridExtra)
library(patchwork)

# Load your data
mutations <- read.csv("FILES/OAPL_n=82_F1CDXRNA-RMC-RET-22-2290_SG44174_21FEB2023145553.txt", sep="\t")
consensus_clusters <- read.csv("RESULTS/Consensus_Clusters.csv")
row.names(consensus_clusters) <- consensus_clusters$X

# Filter mutations to include only those in consensus clusters
mutations <- mutations[mutations$SAMPLE.ID %in% row.names(consensus_clusters),]

# Filter to include only specific variant types
mutations <- mutations[mutations$VARIANT.TYPE %in% c("short-variant", "copy-number-alteration", "rearrangement"),]

# Merge data with consensus clusters and rename column
merged_data <- mutations %>%
  left_join(consensus_clusters, by = c("SAMPLE.ID" = "X")) %>%
  rename(TC_Cluster = final)

# Recode variant type based on CNA.TYPE
merged_data$VARIANT.TYPE <- ifelse(merged_data$CNA.TYPE == "amplification", "copy-number-gain",
                                   ifelse(merged_data$CNA.TYPE == "loss", "copy-number-loss",
                                          merged_data$VARIANT.TYPE))

# Select only relevant columns
cols <- c("SAMPLE.ID", "GENE", "VARIANT.TYPE", "TC_Cluster")
merged_data <- merged_data[, colnames(merged_data) %in% cols]

# Split data by TC_Cluster
split_data <- split(merged_data, merged_data$TC_Cluster)

# Function to create individual plots
create_plot <- function(data, custom_colors) {
  # No need to melt the data frame, directly use the data frame for plotting
  gene_counts <- data %>%
    group_by(GENE, VARIANT.TYPE) %>%
    summarise(Frequency = n()) %>%
    ungroup()
  
  # Create custom order based on frequency
  custom_order <- gene_counts %>%
    group_by(GENE) %>%
    summarise(Total_Frequency = sum(Frequency)) %>%
    filter(Total_Frequency >= 3) %>%
    arrange(desc(Total_Frequency)) %>%
    pull(GENE)
  
  gene_counts$GENE <- factor(gene_counts$GENE, levels = custom_order)
  gene_counts <- gene_counts[!is.na(gene_counts$GENE),]
  
  # Create ggplot with custom order, filtered data, and custom colors
  p <- ggplot(gene_counts, aes(x = GENE, y = Frequency, fill = VARIANT.TYPE)) +
    geom_bar(stat = 'identity', width = 0.8) +  # Adjust bar width
    scale_fill_manual(values = custom_colors) +
    labs(fill = "Variant Type") +  # Change legend title here
    theme_minimal(base_size = 14) +  # Use a minimal theme with adjusted base size
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major.x = element_blank(),  # Remove major gridlines on x-axis
      panel.grid.minor = element_blank(),    # Remove minor gridlines
      panel.grid.major.y = element_line(color = "gray80"),  # Subtle gridlines on y-axis
      plot.margin = margin(10, 10, 10, 10)   # Add margins for better spacing
    ) +
    xlab("Genes") +
    ylab('Frequency') +
    ggtitle(paste("Cluster", data$TC_Cluster[1], "Detected Alterations"))
  
  return(p)
}

# Define custom colors using a more subtle palette
custom_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

# Apply function to each cluster data
plots <- lapply(split_data, create_plot, custom_colors = custom_colors)

# Arrange plots in 2x2 grid and add a single legend at the bottom
combined_plot <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')

# Print the combined plot
print(combined_plot)




