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


clinical_data <- read.csv("RESULTS/TCGA-TC-Classified.csv", row.names = 1)

cinsarc <- read.csv("RESULTS/clinical_cinsarc_paper.csv")

sarculator <- read.csv("RESULTS/Sarculator_TCGA.csv", sep=";")


sarculator$Patient <- gsub("-", ".", sarculator$Patient, fixed = TRUE)


row.names(sarculator) <- sarculator$Patient 

sarculator <- sarculator[row.names(clinical_data),]

clinical_data$SARCULATOR <- sarculator$Sarculator..5.Year.OS.

clinical_data$Age <- sarculator$age

clinical_data$SARCULATOR_STRAT <- ifelse(clinical_data$SARCULATOR <= 0.60,"Low", "High")


extracted_part <- sub("^(([^.]+\\.){3}[^.]{3}).*", "\\1", row.names(clinical_data))
extracted_part <- sub("A$", "", extracted_part)
extracted_part <- sub("B$", "", extracted_part)
new_row_names <- gsub("\\.", "-", extracted_part)

row.names(clinical_data) <- new_row_names

names <- intersect(row.names(clinical_data), cinsarc$X)

row.names(cinsarc) <- cinsarc$X

clinical_data <- clinical_data[names,]
cinsarc <- cinsarc[names,]

clinical_data <- clinical_data[!is.na(clinical_data$Patient),]

clinical_data$CINSARC <- cinsarc$CINSARC



cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~ Histology + PaperHistology + CINSARC  + Cluster + FNCLCC_GRADE,  data = clinical_data)

temp <- cox.zph(cox) 
print(temp)                  # display the results 
plot(temp) 
ggcoxzph(temp)

ggforest(cox, data = clinical_data, fontsize = 1)



all(rownames(clinical_data) == rownames(cinsarc))


clinical_data$Cluster

cinsarc$TC <- clinical_data$Cluster

# Create a survival object using the META_DATE (time) and META (event indicator)
surv_object <- Surv(time = cinsarc$META_DATE, event = cinsarc$META)

# Fit the Kaplan-Meier survival model, stratified by CINSARC group (C1 vs C2)
fit <- survfit(surv_object ~ CINSARC, data = cinsarc)

ggsurvplot(fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis",
           ylab = "Metastatic Free Survival Probability")



# Curves combining C1 and C3, not c1 and not c3, c1 and c3
km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ CINSARC, data=clinical_data)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis")


# Curves combining C1 and C3, not c1 and not c3, c1 and c3
km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Cluster, data=clinical_data)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis")


# Curves combining C1 and C3, not c1 and not c3, c1 and c3
km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Histology, data=clinical_data)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis")

clinical_data$Cluster <- toupper(clinical_data$Cluster)


cm <- table(clinical_data$Cluster, clinical_data$CINSARC)

cm

# Convert the table to a data frame
cm_df <- as.data.frame(cm)

# Plot confusion matrix
ggplot(data = cm_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  labs(x = "CINSARC", y = "Transcriptomic Cluster", fill = "Frequency", title = "Confusion Matrix") +
  theme_minimal()

clinical_data$SARCULATOR <- clinical_data$SARCULATOR_STRAT

# Calculate Cox models and concordance values
cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  CINSARC, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
cinsarc <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  Cluster, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
cluster <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
sarculator <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + CINSARC, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
cinsarculator <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + Cluster, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
tc_sarculator <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + Cluster + CINSARC, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
tc_cinsarculator <- cox$concordance["concordance"]

# New models including age
cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + Age, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
sarculator_age <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  CINSARC + Age, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
cinsarc_age <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  Cluster + Age, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
cluster_age <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + CINSARC + Age, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
cinsarculator_age <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + Cluster + Age, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
tc_sarculator_age <- cox$concordance["concordance"]

cox <- coxph(Surv(Last_FU, as.numeric(Status)) ~  SARCULATOR + Cluster + CINSARC + Age, data = clinical_data)
#ggforest(cox, data = clinical_data, fontsize = 1)
tc_cinsarculator_age <- cox$concordance["concordance"]

# Compile all concordance values into a data frame
concordance_values <- data.frame(
  Model = c("SARCULATOR", "SARCULATOR + CINSARC", "SARCULATOR + TC", 
            "SARCULATOR + TC + CINSARC", "CINSARC", "TC",
            "SARCULATOR + Age", "CINSARC + Age", "TC + Age", 
            "SARCULATOR + CINSARC + Age", "SARCULATOR + TC + Age", 
            "SARCULATOR + TC + CINSARC + Age"),
  Concordance = c(sarculator, cinsarculator, tc_sarculator, 
                  tc_cinsarculator, cinsarc, cluster,
                  sarculator_age, cinsarc_age, cluster_age,
                  cinsarculator_age, tc_sarculator_age, tc_cinsarculator_age)
)

# Order the data frame by Concordance values
concordance_values <- concordance_values[order(concordance_values$Concordance, decreasing = TRUE), ]


concordance_values$Color <- ifelse(grepl("TC", concordance_values$Model), "orange", "skyblue")

concordance_values$Label <- ifelse(grepl("TC", concordance_values$Model), "Includes TC", "Does Not Include TC")

# Plot the ordered barplot with conditional coloring

# Plot the ordered barplot with conditional coloring and a legend
ggplot(concordance_values, aes(x = reorder(Model, Concordance), y = Concordance, fill = Label)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Concordance, 2)), hjust = -0.1, size = 4) +
  coord_flip() +
  xlab("Model") +
  ylab("Concordance") +
  ggtitle("Concordance Indexes of Different Models") +
  theme_minimal() +
  scale_fill_manual(values = c("Includes TC" = "darkorange1", "Does Not Include TC" = "skyblue"))  

table(clinical_data$FNCLCC_GRADE)

clinical_data_stra_1 <- clinical_data[clinical_data$SARCULATOR == "Low",]
clinical_data_stra_2 <- clinical_data[clinical_data$SARCULATOR == "High",]

km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Cluster, data=clinical_data_stra_1)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis",
           title = "Overall Survival Sarculator Score <= 0.60")

km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Cluster, data=clinical_data_stra_2)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis",
           title = "Overall Survival Sarculator Score > 0.60")

clinical_data$strat <- ifelse(clinical_data$SARCULATOR <= 0.60, "\"<=0.60\" ", "\">0.60\"")

km_fit <- survfit(Surv(Last_FU, as.numeric(Status)) ~ Cluster, data=clinical_data)

ggsurvplot(km_fit,pval=TRUE,
           risk.table=TRUE,
           conf.int = TRUE, 
           legend.title = "Survival Analysis")
