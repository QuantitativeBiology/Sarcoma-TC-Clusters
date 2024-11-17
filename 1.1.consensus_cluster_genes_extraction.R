# Load required libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(edgeR)
library(dplyr)
library(readxl)
library(EnhancedVolcano)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(tibble)
library(gridExtra)
library(ComplexHeatmap)
library(VennDiagram)

DNA_data <- read_excel("FILES/OSPL_n=82_F1CDXRNA-RMC-RET-22-2290_SG44174_21FEB2023145617.xlsx")

RNA_data <- read.csv("FILES/22-2290 Roche RNA GEP-COUNTS 01MAY2023.txt", sep = "\t")

row.names(RNA_data) = sub("\\*.*", "", RNA_data$SPECIMEN)

cluster_data <- read.csv("RESULTS/Consensus_Clusters.csv",row.names =1)

cluster_data$extra = cluster_data$final

cluster_data <- cluster_data[order(cluster_data$final),]

# Remove unnecesssary collumns
RNA_data <- subset(RNA_data, select = -c(QC_STATUS, BAITSET, TEST_TYPE, EXPRESSION_UNIT, QC_FLAGS))

# Remove last row # FOR RESEARCH USE ONLY
RNA_data <- RNA_data[-nrow(RNA_data), ]

Fmi_sample_id <- DNA_data[,c("FMI SAMPLE ID", "SAMPLE ID")]

Fmi_sample_id <- DNA_data[Fmi_sample_id$`FMI SAMPLE ID` %in% row.names(RNA_data),]

order_indices <- match(Fmi_sample_id$`FMI SAMPLE ID`,row.names(RNA_data))

RNA_data <- RNA_data[order_indices,]

all(Fmi_sample_id$`FMI SAMPLE ID` == row.names(RNA_data))

row.names(RNA_data) = Fmi_sample_id$`SAMPLE ID`

all(row.names(RNA_data) == Fmi_sample_id$`SAMPLE ID`)

RNA_data <- RNA_data[row.names(RNA_data) %in% row.names(cluster_data),]
cluster_data <- cluster_data[rownames(cluster_data) %in% row.names(RNA_data),,drop=FALSE]

all(row.names(cluster_data) == row.names(RNA_data))

order_indices <- match(row.names(cluster_data),row.names(RNA_data))
RNA_data <- RNA_data[order_indices,]

all(row.names(cluster_data) == row.names(RNA_data))

RNA_data$SPECIMEN = NULL
RNA_data <- t(RNA_data)
RNA_data <- RNA_data[rowSums(RNA_data) > 0, ]

RNA_data <- RNA_data[,colnames(RNA_data) %in% row.names(cluster_data) ]

# make code above this line more simple


cluster_one <- as.numeric(grepl("1", cluster_data$final, ignore.case = TRUE))

cluster_two <- as.numeric(grepl("2", cluster_data$final, ignore.case = TRUE))

cluster_three <- as.numeric(grepl("3", cluster_data$final, ignore.case = TRUE))

cluster_four <- as.numeric(grepl("4", cluster_data$final, ignore.case = TRUE))


RNA_data <- RNA_data[,row.names(cluster_data)]

all(row.names(cluster_data) == colnames(RNA_data))

design <- cbind(cluster_one,cluster_two, cluster_three, cluster_four)

keep <- filterByExpr(RNA_data, design = design)

RNA_data <- RNA_data[keep,]

RNA_data <- DGEList(counts = RNA_data, genes = rownames(RNA_data))

cont.matrix <- makeContrasts(cluster_one - cluster_two ,
                             cluster_one - cluster_three,
                             cluster_one - cluster_four, 
                             cluster_two - cluster_three,
                             cluster_two - cluster_four,
                             cluster_three - cluster_four,
                             levels=design)

Voom <- voom(RNA_data, design, plot = FALSE, normalize.method = "quantile") 

vfit <- lmFit(Voom, design) 

vfit  <- contrasts.fit(vfit,cont.matrix) 

efit <- eBayes(vfit) 

summary(decideTests(efit))
deg <- topTable(efit, coef =1, adjust.method = 'fdr', number = Inf)
deg2 <- topTable(efit, coef =2, adjust.method = 'fdr', number = Inf)
deg3 <- topTable(efit, coef =3, adjust.method = 'fdr', number = Inf)
deg4 <- topTable(efit, coef =4, adjust.method = 'fdr', number = Inf)
deg5 <- topTable(efit, coef =5, adjust.method = 'fdr', number = Inf)
deg6 <- topTable(efit, coef =6, adjust.method = 'fdr', number = Inf)

x <- rownames(deg)
y <- deg$adj.P.Val


for (i in 1:length(rownames(deg))) {
  if (y[i] >= 0.05) {
    x[i] <- ""
  }
}

x1 <- EnhancedVolcano(
  deg,
  xlim = c(-8,8),
  ylim = c(0,30),
  lab = x,
  x = 'logFC',
  y = 'P.Value',
  labSize = 3.0,
  pCutoff = 1e-02,
  FCcutoff = 0.75,
  title = 'Cluster One vs Cluster Two',
  legendPosition = 'none',
  drawConnectors = FALSE,
  max.overlaps = 100
)
#
#
x <- rownames(deg2)
y <- deg2$adj.P.Val

for (i in 1:length(rownames(deg2))) {
  if (y[i] >= 0.05) {
    x[i] <- ""
  }
}

x2 <- EnhancedVolcano(
  deg2,
  lab = x,
  xlim = c(-8,8),
  ylim = c(0,30),
  x = 'logFC',
  y = 'P.Value',
  labSize = 3.0,
  pCutoff = 1e-02,
  FCcutoff = 0.75,
  title = 'Cluster One vs Cluster Three',
  legendPosition = 'none',
  drawConnectors = FALSE,
  max.overlaps = 100
)
#
x <- rownames(deg3)
y <- deg3$adj.P.Val
#
for (i in 1:length(rownames(deg3))) {
  if (y[i] >= 0.05) {
    x[i] <- ""
  }
}

x3 <- EnhancedVolcano(
  deg3,
  lab = x,
  xlim = c(-8,8),
  ylim = c(0,30),
  x = 'logFC',
  y = 'P.Value',
  labSize = 3.0,
  pCutoff = 1e-02,
  FCcutoff = 0.75,
  title = 'Cluster One vs Cluster Four',
  legendPosition = 'none',
  drawConnectors = FALSE,
  max.overlaps = 100
)
#
x <- rownames(deg4)
y <- deg4$adj.P.Val

for (i in 1:length(rownames(deg4))) {
  if (y[i] >= 0.05) {
    x[i] <- ""
  }
}
#
x4 <- EnhancedVolcano(
  deg4,
  lab = x,
  xlim = c(-8,8),
  ylim = c(0,30),
  x = 'logFC',
  y = 'P.Value',
  labSize = 3.0,
  pCutoff = 1e-02,
  FCcutoff = 0.75,
  title = 'Cluster Two vs Cluster Three',
  legendPosition = 'none',
  drawConnectors = FALSE,
  max.overlaps = 100
)
#
#
x <- rownames(deg5)
y <- deg5$adj.P.Val

for (i in 1:length(rownames(deg5))) {
  if (y[i] >= 0.05) {
    x[i] <- ""
  }
}

x5 <- EnhancedVolcano(
  deg5,
  lab = x,
  xlim = c(-8,8),
  ylim = c(0,30),
  x = 'logFC',
  y = 'P.Value',
  labSize = 3.0,
  pCutoff = 1e-02,
  FCcutoff = 0.75,
  title = 'Cluster Two vs Cluster Four',
  legendPosition = 'none',
  
  drawConnectors = FALSE,
  max.overlaps = 100
)
#
#
x <- rownames(deg6)
y <- deg6$adj.P.Val

for (i in 1:length(rownames(deg6))) {
  if (y[i] >= 0.05) {
    x[i] <- ""
  }
}
#
x6 <- EnhancedVolcano(
  deg6,
  lab = x,
  xlim = c(-8,8),
  ylim = c(0,30),
  x = 'logFC',
  y = 'P.Value',
  labSize = 3.0,
  pCutoff = 1e-02,
  FCcutoff = 0.75,
  title = 'Cluster Three vs Cluster Four',
  legendPosition = 'none',
  drawConnectors = FALSE,
  max.overlaps = 100
)


gridExtra::grid.arrange(x1,x2,x3, ncol=3)
gridExtra::grid.arrange(x4,x5,x6, ncol=3)


deg <- deg[deg$adj.P.Val< 0.05 ,]
deg2 <- deg2[deg2$adj.P.Val< 0.05 ,]
deg3 <- deg3[deg3$adj.P.Val< 0.05 ,]
deg4 <- deg4[deg4$adj.P.Val< 0.05 ,]
deg5 <- deg5[deg5$adj.P.Val< 0.05 ,]
deg6 <- deg6[deg6$adj.P.Val< 0.05 ,]

c1 <- intersect(deg$genes,deg2$genes)
c1 <- intersect(c1, deg3$genes)

c2 <- intersect(deg$genes, deg4$genes)
c2 <- intersect(c2,deg5$genes)

c3 <- intersect(deg2$genes, deg4$genes)
c3 <- intersect(c3,deg6$genes)

c4 <- intersect(deg3$genes,deg5$genes)
c4 <- intersect(c4,deg6$genes)

all_genes <- unique(c(c1, c2, c3, c4))


# Get over expressed on cluster 1
deg_c1 <- deg[deg$logFC > 0.00,]
deg_c1_1 <- deg2[deg2$logFC > 0.00,]
deg_c1_2 <- deg3[deg3$logFC > 0.00,]

c1_over <- intersect(deg_c1$genes,deg_c1_1$genes)
c1_over <- intersect(c1_over, deg_c1_2$genes)

#under

deg_c1 <- deg[deg$logFC < 0.00,]
deg_c1_1 <- deg2[deg2$logFC < 0.00,]
deg_c1_2 <- deg3[deg3$logFC < 0.00,]

c1_under <- intersect(deg_c1$genes,deg_c1_1$genes)
c1_under <- intersect(c1_under, deg_c1_2$genes)


# Get over expressed on cluster 2
deg_c2 <- deg[deg$logFC < 0.00,]
deg_c2_1 <- deg4[deg4$logFC > 0.00,]
deg_c2_2 <- deg5[deg5$logFC > 0.00,]

c2_over <- intersect(deg_c2$genes,deg_c2_1$genes)
c2_over <- intersect(c2_over, deg_c2_2$genes)

#under

deg_c2 <- deg[deg$logFC > 0.00,]
deg_c2_1 <- deg4[deg4$logFC < 0.00,]
deg_c2_2 <- deg5[deg5$logFC < 0.00,]

c2_under <- intersect(deg_c2$genes,deg_c2_1$genes)
c2_under <- intersect(c2_under, deg_c2_2$genes)



# Get over expressed on cluster 3
deg_c3 <- deg2[deg2$logFC < 0.00,]
deg_c3_1 <- deg4[deg4$logFC < 0.00,]
deg_c3_2 <- deg6[deg6$logFC > 0.00,]

c3_over <- intersect(deg_c3$genes,deg_c3_1$genes)
c3_over <- intersect(c3_over, deg_c3_2$genes)

#under

deg_c3 <- deg2[deg2$logFC > 0.00,]
deg_c3_1 <- deg4[deg4$logFC > 0.00,]
deg_c3_2 <- deg6[deg6$logFC < 0.00,]

c3_under <- intersect(deg_c3$genes,deg_c3_1$genes)
c3_under <- intersect(c3_under, deg_c3_2$genes)


# Get over expressed on cluster 4
deg_c4 <- deg3[deg3$logFC < 0.00,]
deg_c4_1 <- deg5[deg5$logFC < 0.00,]
deg_c4_2 <- deg6[deg6$logFC < 0.00,]

c4_over <- intersect(deg_c4$genes,deg_c4_1$genes)
c4_over <- intersect(c4_over, deg_c4_2$genes)

#under

deg_c4 <- deg3[deg3$logFC > 0.00,]
deg_c4_1 <- deg5[deg5$logFC > 0.00,]
deg_c4_2 <- deg6[deg6$logFC > 0.00,]

c4_under <- intersect(deg_c4$genes,deg_c4_1$genes)
c4_under <- intersect(c4_under, deg_c4_2$genes)


list_names <- c("c1_over", "c1_under", "c2_over", "c2_under", "c3_over", "c3_under", "c4_over", "c4_under")

# Create a vector with the corresponding lists
lists <- list(c1_over, c1_under, c2_over, c2_under, c3_over, c3_under, c4_over, c4_under)

# Create a file connection for writing
file_conn <- file("RESULTS/Consensus_Clusters_Genes.txt", open = "wt")

# Loop through the lists and write them to the file
for(i in seq_along(lists)) {
  cat(paste(list_names[i], ":\n", toString(lists[[i]]), "\n\n"), file = file_conn)
}

# Close the file connection
close(file_conn)


















