library(ConsensusClusterPlus)
library(readxl)
library(ComplexHeatmap)
library(dplyr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(GGally)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

tc_clusters <- read.csv("RESULTS/Consensus_Clusters.csv")

# File path to the gene lists
file_path <- "RESULTS/Consensus_Clusters_Genes.txt"

row.names(tc_clusters) <- tc_clusters$X

tc_clusters$final <- ifelse(tc_clusters$final == 1, "C1" ,tc_clusters$final)
tc_clusters$final <- ifelse(tc_clusters$final == 2, "C2" ,tc_clusters$final)
tc_clusters$final <- ifelse(tc_clusters$final == 3, "C3" ,tc_clusters$final)
tc_clusters$final <- ifelse(tc_clusters$final == 4, "C4" ,tc_clusters$final)

clinical_data_clean <- read.csv("FILES/clinical_data_clean.csv", row.names = 1)

voom <- read.csv("FILES/normalized_counts.csv", row.names = 1)
colnames(voom) <- sub("^X","", colnames(voom))

common <- intersect(colnames(voom), row.names(clinical_data_clean))

tc_clusters <- tc_clusters[common,]
tc_clusters <- tc_clusters[order(tc_clusters$final),]

clinical_data_clean <- clinical_data_clean[row.names(tc_clusters),]
voom <- voom[,row.names(tc_clusters)]

read_gene_lists <- function(file_path) {
  # Read all lines from the file
  lines <- readLines(file_path)
  
  # Initialize variables
  lists <- list()
  current_list_name <- NULL
  
  # Iterate through each line
  for (line in lines) {
    line <- trimws(line)  # Remove extra whitespace from the line
    
    if (grepl(":", line)) {  # If the line defines a list name
      current_list_name <- sub(":", "", line)  # Extract list name
      current_list_name <- trimws(current_list_name)  # Remove any extra spaces from the name
      lists[[current_list_name]] <- character(0)  # Initialize the list
    } else if (nchar(line) > 0 && !is.null(current_list_name)) {
      # Append elements to the current list
      lists[[current_list_name]] <- c(lists[[current_list_name]], unlist(strsplit(line, ",\\s*")))
    }
  }
  
  # Assign each list to a variable in the global environment
  for (name in names(lists)) {
    assign(name, lists[[name]], envir = .GlobalEnv)
  }
  
  return(lists)
}

# Call the function to read and assign variables
gene_lists <- read_gene_lists(file_path)

united_gene_sets <- list(
  c1_over = c1_over,
  c1_under = c1_under,
  c2_over = c2_over,
  c2_under = c2_under,
  c3_over = c3_over,
  c3_under = c3_under,
  c4_over = c4_over,
  c4_under = c4_under
)


# Find the common genes between the united_gene_sets and the row names of voom
common_genes <- intersect(unlist(united_gene_sets), rownames(voom))


# Filter each gene list to keep only the common genes
filtered_gene_sets <- lapply(united_gene_sets, function(genes) intersect(genes, common_genes))

# Combine the filtered gene sets into a single vector
filtered_gene_list <- unlist(filtered_gene_sets)

# Generate the vector with the names of the gene sets repeated according to their lengths
gene_set_names <- unlist(lapply(names(filtered_gene_sets), function(name) rep(name, length(filtered_gene_sets[[name]]))))

# Ensure voom only contains the filtered genes
voom_filtered <- voom[filtered_gene_list, ]

voom_original <- voom

voom <- voom_filtered

voom <- data.frame(voom)

voom <- voom[unlist(united_gene_sets), ]

gene_set_colors <- c("c1_over" = "#4B0082",  # Dark purple
                     "c1_under" = "#9370DB",  # Light purple
                     "c2_over" = "#2F4F4F",  # Dark grey
                     "c2_under" = "#D3D3D3",  # Light grey
                     "c3_over" = "#FF1493",  # Dark pink
                     "c3_under" = "#FFB6C1",  # Light pink
                     "c4_over" = "#FFD700",  # Dark yellow
                     "c4_under" = "#FFFFE0")  # Light yellow

# Create a row annotation object for the gene clusters
row_ha <- rowAnnotation(
  GeneSet = factor(gene_set_names, levels = unique(gene_set_names)),
  col = list(GeneSet = gene_set_colors)
)

column_ha = HeatmapAnnotation(
  Transcriptomic = tc_clusters$final,
  Histopathological = clinical_data_clean$Sarcoma.Histopathological.Subtype,
  Gender = clinical_data_clean$Sexo,
  Age = clinical_data_clean$Age,
  Local_Recurrence = clinical_data_clean$Recidiva.S.N,
  Distant_Recurrence = clinical_data_clean$Metastização.S.N,
  col = list(Transcriptomic = c("C1" = "purple",
                                "C2" = "grey",
                                "C3" = "pink",
                                "C4" = "yellow"),
             Histopathological = c("DDLPS" = "red",
                                          "LMS" = "green",
                                          "UPS" = "lightblue"),
  Local_Recurrence = c("1"= "black",
                       "0" = "white"),
  Distant_Recurrence = c("1"= "black",
                         "0" = "white")
  )
)


voom <- data.frame(t(voom))

#voom_needed <- data.frame(t(voom_needed))

#voom$CDK4 <- voom_needed$CDK4
#voom$CDKN2A <- voom_needed$CDKN2A
#voom$BRCA2 <- voom_needed$BRCA2


voom <- scale(voom)

voom <- data.frame(t(voom))


Heatmap(as.matrix(voom), 
        name = "Expression", 
        top_annotation = column_ha,
        right_annotation = row_ha,
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        row_title = "Normalized Gene Expression", 
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 5),  # Adjust the font size
        column_names_gp = gpar(fontsize = 10),
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE)


voom_original <- data.frame(t(voom_original))

voom_original$Cluster <- tc_clusters$final

# Extract the relevant columns
cdk4_data <- voom_original[, c("HLA.DMA", "Cluster")]

# Create the violin plot
ggplot(cdk4_data, aes(x = factor(Cluster), y = HLA.DMA)) +
  geom_violin(trim = FALSE, fill = "lightgrey", color = "black") +  # Violin plot
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +  # Overlay boxplot
  geom_jitter(width = 0.2, size = 1, color = "black", alpha = 0.7) +  # Add data points
  labs(title = "HLA.DMA Expression Across Transcriptomic Clusters",
       x = "Transcriptomic Cluster",
       y = "HLA.DMA Expression (Voom Normalized)") +
  theme_classic(base_size = 14) +  # Classic theme with larger font
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered bold title
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(color = "black"),  # Black axis text
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Border around plot
  )


