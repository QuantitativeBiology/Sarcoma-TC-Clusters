library(EnhancedVolcano)
library(limma)
library(ComplexHeatmap)
library(gplots)
library(edgeR)
library(fgsea)
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

c1_over <- c("CCND2", "SOX18", "CBFA2T3", "PAX5", "CCL19", "CX3CL1", "PNOC", "OLFM1", 
             "BCAM", "ECSCR", "PLVAP", "CD79B", "TTYH1", "GRM4", "CDK4", "ERG", "DTX1", "FMOD", "ZBTB46", "NKD1")

c1_under <- c("RRAGC", "RAC1", "CCNB2", "KIF23", "NEK2", "PBK", "GMNN", "AURKA", "CCNA2", "WHSC1", "NUF2", "CDCA8", "UBE2T", "XPO1", "ANLN", "ECT2", "MALT1", "RRM2", "CRNDE", "CCNB1", "CDKN3", "CDC25C", "TOP2A", "STIL", "BUB1B", "CHEK1", "TTK", "MCM4", "CDCA5", "BUB1", "EXO1", "CEP55", "CDC20", "UBE2C", "ABL2", "DEK", "BRCA1", "FGFR1OP", "PTTG1", "SNW1", "PALB2", "CENPF", "RAD51AP1", "AURKB", "NDC80", "TPX2", "TYMS", "RAD54L", "FANCD2", "BRIP1", "MELK", "GINS2", "CENPM", "KIF2C", "RAD51", "HIST1H3B", "BRCA2", "EPS15")

c2_over <- c("MAGEC2", "SSX2", "SSX1", "MAGEA3", "MAGEA12", "SSX3", "MAGEB2", "SSX2B", "MAGEA2B", "FZD6", "MAGEB1", "JAZF1", "ACVR1C", "BAP1", "CTNNB1", "RGS16", "MMP11", "MRAS")

c2_under <- c("DHX58", "IL12A", "TNFRSF1B", "MAP3K8", "SULT1A1", "LILRB5")

c3_over <- c("HAVCR2", "FCGR3B", "CYBB", "LCP1", "CCL2", "HLA.DMA", "FPR3", "ITGB2", "CXCL10", "CD84", "HLA.DMB", "CSF2", "HMGA1", "BTK", "HLA.DRA", "IL21R", "GMFG", "CYLD", "ATIC", "CCL18", "LAIR1", "PLEK2", "CCR5", "JAML", "IL7R", "RGS10", "HLA.DRB1", "SYK", "FAM26F", "BATF", "NFKB2", "CD74", "PSMB10", "ETV5", "TGFB1", "FGR", "SERPINE1", "SEMA7A", "CD3G", "HLA.DQA1", "GBP5", "PHF11", "HLA.DOA", "FN1", "KCNMA1")

c3_under <- c("TCF7L1", "PBX1", "SMAD9", "FZD7", "BCL9", "TP53INP2", "SCUBE2", "NRTN", "DHH", "TET1", "FOXO1", "FGFR3", "HAP1", "AXIN2", "FOXO6", "HES1", "BMP4", "GAS1", "CDKN1C", "SESN3", "PDGFD", "LINC00598", "ZNF521", "SEMA6D", "CITED4", "SH3PXD2A", "WNT11", "GPC4", "TMEM38A", "TRIM2", "DOT1L", "FOXC1", "DCLK1", "PRKACG", "KDM5C", "FGFR2", "PHLPP1", "FOXO4", "TCF7L2", "CDKN2A")

c4_over <- c("PHLPP1", "IL6ST", "DDR2", "TP73", "FAM64A", "SYCP3", "CLDN4", "LINC.ROR", "CLCA2", "SMAD3", "GAS7", "CFD", "PDGFD", "ADRB2", "PGR", "CD34", "TEK")

c4_under <- c("HOPX", "FZD2", "NBEAP1", "ACTN1")

cinsarc_genes <-  c("MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2", "ESR1", "PGR", "BCL2", 
                    "SCUBE2", "MMP11", "CTSL2", "GRB7", "ERBB2", "ACTB", "GAPDH", "GUS", 
                    "RPLP0", "TFRC", "GSTM1", "CD68", "BAG1")

united_gene_sets <- list(
  c1 = c(c1_over, c1_under),
  c2 = c(c2_over, c2_under),
  c3 = c(c3_over, c3_under),
  c4 = c(c4_over, c4_under),
  cinsarc = cinsarc_genes
)

gene_sets <- list(
  c1_over = c1_over,
  c1_under = c1_under,
  c2_over = c2_over,
  c2_under = c2_under,
  c3_over = c3_over,
  c3_under = c3_under,
  c4_over = c4_over,
  c4_under = c4_under,
  CINSARC = cinsarc_genes
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

# Shime1x_syn_spearman1 <- Heatmap(as.matrix(nes_mat),
#                                  heatmap_legend_param=list(title="NES", 
#                                                            direction='vertical'),
#                                  cluster_rows = TRUE,
#                                  cluster_columns = TRUE,
#                                  column_names_gp = gpar(col = NA)
# 
#                                  )
# sx_syn1_spear <- draw(Shime1x_syn_spearman1, heatmap_legend_side='left')



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
#
clinical_data <- clinical_data[clinical_data$Cluster != "c2",]
clinical_data <- clinical_data[clinical_data$Cluster != "c4",]

adjusted_p_values_mat <- adjusted_p_values_mat[row.names(adjusted_p_values_mat) %in% row.names(clinical_data),]


clinical_data <- clinical_data[row.names(clinical_data) %in% row.names(adjusted_p_values_mat),]




cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~Histology + Cluster + PaperHistology + FNCLCC_GRADE,  data = clinical_data)


ggforest(cox,data=clinical_data, fontsize = 1)


anova(cox)


km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Cluster, data=clinical_data)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis")

write.csv(clinical_data,"RESULTS/TCGA-TC-Classified.csv")





