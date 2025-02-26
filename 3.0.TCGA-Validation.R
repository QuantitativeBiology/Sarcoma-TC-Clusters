library(EnhancedVolcano)
library(limma)
library(ComplexHeatmap)
library(gplots)
library(edgeR)
library(corto)
library(survival)
library(survminer)
library("PerformanceAnalytics")
library(msigdbr)
library(dplyr)

RNA_data <- read.csv("RESULTS/TCGA-SARC.csv",row.names=1)

clinical_data <- read.csv("RESULTS/TCGA-SARC-CLINICAL.csv",row.names=1)

cinsarc <- read.csv("RESULTS/clinical_cinsarc_paper.csv")

sarculator <- read.csv("RESULTS/Sarculator_TCGA.csv", sep=";")

file_path <- "RESULTS/Consensus_Clusters_Genes.txt"

clinical_data <- clinical_data[clinical_data$Histology %in% c("Leiomyosarcoma, NOS", "Dedifferentiated liposarcoma"),]

clinical_data$Patient <- gsub("-", ".", clinical_data$Patient, fixed = TRUE)

clinical_data <- clinical_data[!is.na(clinical_data$Patient),]

clinical_data$FNCLCC_GRADE <- as.character(clinical_data$FNCLCC_GRADE)

RNA_data <- RNA_data[, colnames(RNA_data) %in% clinical_data$Patient]

clinical_data <- clinical_data[order(clinical_data$PaperHistology),]

RNA_data <- RNA_data[,clinical_data$Patient]

keep <- filterByExpr(RNA_data)

RNA_data <- RNA_data[keep,]

Voom <- voom(RNA_data, plot = FALSE,normalize.method = "quantile")

normalized <- Voom$E

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
  c1 = c(c1_over, c1_under),
  c2 = c(c2_over, c2_under),
  c3 = c(c3_over, c3_under),
  c4 = c(c4_over, c4_under)
)

gene_sets <- list(
  c1_over = c1_over,
  c1_under = c1_under,
  c2_over = c2_over,
  c2_under = c2_under,
  c3_over = c3_over,
  c3_under = c3_under,
  c4_over = c4_over,
  c4_under = c4_under
)

# single sample GSEA

hs_kegg_df <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(
    gs_cat == "C2",
    gs_subcat == "CP:KEGG"
  ) %>%
  dplyr::select(gs_name, gene_symbol)

# Create a list of pathways and associated genes
gene_sets_all <- hs_kegg_df %>%
  group_by(gs_name) %>%
  summarize(genes = list(unique(gene_symbol))) %>%
  ungroup()

# Convert the result to a named list
gene_sets_all <- setNames(gene_sets_all$genes, gene_sets_all$gs_name)

combined_gene_sets <- c(gene_sets_all, gene_sets)

# Apply ssGSEA over all genesets
res <- corto::ssgsea(normalized,combined_gene_sets,minsize = 3)

nes_mat <- as.data.frame(res)

# Obtain p values from NES
p_val_res <- corto::z2p(res)


adjusted_p_values_mat <- matrix(NA, nrow = nrow(p_val_res), ncol = ncol(p_val_res))

# Adjust p values 
adjusted_p_values_mat <- p.adjust(as.vector(p_val_res), method = "BH")
adjusted_p_values_mat <- matrix(adjusted_p_values_mat, nrow = nrow(p_val_res), ncol = ncol(p_val_res))

adjusted_p_values_mat <- data.frame(adjusted_p_values_mat)

# Keep only out Gene Sets
colnames(p_val_res) <- colnames(normalized)

row.names(adjusted_p_values_mat) <- row.names(nes_mat)

adjusted_p_values_mat <- adjusted_p_values_mat[names(gene_sets),]

colnames(adjusted_p_values_mat) <- colnames(normalized)


nes_mat <- nes_mat[names(gene_sets),]

# Observe Correlation between NES

#chart.Correlation.nostars(t(nes_mat), histogram=TRUE, method="spearman")

# adjusted_p_values_mat , nes_mat, plot heatmap with p values

 Shime1x_syn_spearman1 <- Heatmap(as.matrix(nes_mat),
                                  heatmap_legend_param=list(title="NES", 
                                                            direction='vertical'),
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  column_names_gp = gpar(col = NA)
 
                                  )
sx_syn1_spear <- draw(Shime1x_syn_spearman1, heatmap_legend_side='left')



cluster_assignments_p_value <- data.frame(patient_id = colnames(adjusted_p_values_mat), cluster = NA)

#### Classification according to min adjusted p value

# Iterate over each patient
for (col in colnames(adjusted_p_values_mat)) {
  # Find the row with the maximum absolute enrichment score
  max_row <- which.min(abs(adjusted_p_values_mat[, col]))
  
  # Assign the cluster to the data frame
  cluster_assignments_p_value[cluster_assignments_p_value$patient_id == col, "cluster"] <- rownames(adjusted_p_values_mat)[max_row]
}

cluster_assignments_p_value$cluster <- gsub("_.*", "", cluster_assignments_p_value$cluster)


row.names(clinical_data) <- clinical_data$Patient

clinical_data <- clinical_data[cluster_assignments_p_value$patient_id,]

clinical_data$Cluster <- cluster_assignments_p_value$cluster

adjusted_p_values_mat <- data.frame(t(adjusted_p_values_mat))


colnames(nes_mat) <- colnames(normalized)

p_val_res <- data.frame(t(p_val_res))

clinical_data <- clinical_data[!is.na(clinical_data$Last_FU),]

clinical_data<-clinical_data[clinical_data$Cluster %in% c("c1","c3"),]

adjusted_p_values_mat <- adjusted_p_values_mat[row.names(adjusted_p_values_mat) %in% row.names(clinical_data),]

clinical_data <- clinical_data[row.names(clinical_data) %in% row.names(adjusted_p_values_mat),]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~Histology +PaperHistology + FNCLCC_GRADE + Cluster,  data = clinical_data)

ggforest(cox,data=clinical_data, fontsize = 1)

anova(cox)

colnames(clinical_data)[colnames(clinical_data) == "Cluster"] <- "Transcriptomic_Cluster"

km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Transcriptomic_Cluster, data=clinical_data)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis")

colnames(clinical_data)[colnames(clinical_data) == "Transcriptomic_Cluster"] <- "Cluster"

write.csv(clinical_data,"RESULTS/TCGA-TC-Classified.csv")
