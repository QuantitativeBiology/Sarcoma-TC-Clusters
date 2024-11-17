library(tidyverse)
library(survival)
library(survminer)
library(readxl)
library(forestploter)


clinical_data_clean <- read.csv("FILES/clinical_data_clean.csv", row.names=1)


colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Sarcoma.Histopathological.Subtype"] <- "Sarcoma Histopathological Subtype"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Neo.Adjuvant...Adjuvant.Treatment"] <- "Neo Adjuvant / Adjuvant Treatment"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Local.Recurrence"] <- "Local Recurrence"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Distant.Recurrence"] <- "Distant Recurrence"


sarculator <- read_excel("RESULTS/Sarculator_Results.xlsx")

consensus_clusters <- read.csv("RESULTS/Consensus_Clusters.csv")
row.names(consensus_clusters) <- consensus_clusters$X

names <- consensus_clusters$X


consensus_clusters <- consensus_clusters[names,]

clinical_data_clean <- clinical_data_clean[row.names(clinical_data_clean) %in% names,]
clinical_data_clean <- clinical_data_clean[names,]


row.names(sarculator) <- sarculator$...1

consensus_clusters <- consensus_clusters[row.names(clinical_data_clean),]

clinical_data_clean$Transcriptomic <- consensus_clusters$final


clinical_data_clean$Transcriptomic <- ifelse(clinical_data_clean$Transcriptomic == 1, "C1" ,clinical_data_clean$Transcriptomic)
clinical_data_clean$Transcriptomic <- ifelse(clinical_data_clean$Transcriptomic == 2, "C2" ,clinical_data_clean$Transcriptomic)
clinical_data_clean$Transcriptomic <- ifelse(clinical_data_clean$Transcriptomic == 3, "C3" ,clinical_data_clean$Transcriptomic)
clinical_data_clean$Transcriptomic <- ifelse(clinical_data_clean$Transcriptomic == 4, "C4" ,clinical_data_clean$Transcriptomic)
clinical_data_clean <- clinical_data_clean[!is.na(clinical_data_clean$Nº.Histologia),]

table <- data.frame(row.names = names, pvalue = numeric(length(names)))


clinical_data_clean$`Neo Adjuvant / Adjuvant Treatment` <- as.character(clinical_data_clean$`Neo Adjuvant / Adjuvant Treatment`)
clinical_data_clean$`Local Recurrence` <- as.character(clinical_data_clean$`Local Recurrence`)

sarculator <- sarculator[row.names(clinical_data_clean),]
clinical_data_clean$Sarculator_5_Year_OS <- as.numeric(sarculator$Sarculator)


colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Transcriptomic"] <- "Transcriptomic Clusters"

clinical_data_clean$Sarculator_5_Year_OS <- ifelse(clinical_data_clean$Sarculator_5_Year_OS <= 0.60, "Low", "High")


cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Gender + `Sarcoma Histopathological Subtype` + `Transcriptomic Clusters` + `Neo Adjuvant / Adjuvant Treatment` + `Local Recurrence`  , data = clinical_data_clean)

ggforest(cox,data = clinical_data_clean, fontsize = 1.5)


temp <- cox.zph(cox) 
print(temp)                  # display the results 
plot(temp) 
ggcoxzph(temp)
anova(cox)

clinical_data_clean <- clinical_data_clean[!is.na(clinical_data_clean$Sarculator_5_Year_OS),]

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  `Transcriptomic Clusters` , data = clinical_data_clean)
ggforest(cox, data = clinical_data_clean, fontsize = 1)
tc  <- cox$concordance["concordance"]

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Sarculator_5_Year_OS, data = clinical_data_clean)
sarculator <- cox$concordance["concordance"]

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Age, data = clinical_data_clean)
Age <- cox$concordance["concordance"]

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Sarculator_5_Year_OS + `Transcriptomic Clusters` , data = clinical_data_clean)
tc_sarculator <- cox$concordance["concordance"]

# New models including age
cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Sarculator_5_Year_OS + Age, data = clinical_data_clean)
sarculator_age <- cox$concordance["concordance"]

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  `Transcriptomic Clusters`  + Age, data = clinical_data_clean)
tc_age <- cox$concordance["concordance"]

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Sarculator_5_Year_OS + `Transcriptomic Clusters`  + Age, data = clinical_data_clean)
tc_sarculator_age <- cox$concordance["concordance"]

# Compile all concordance values into a data frame
concordance_values <- data.frame(
  Model = c("SARCULATOR", "TC + SARCULATOR", 
            "TC",
            "Age",
            "TC + Age"),
  Concordance = c(sarculator, tc_sarculator, 
                  tc, Age,
                  tc_age)
)

concordance_values <- concordance_values[order(concordance_values$Concordance, decreasing = TRUE), ]


concordance_values$Color <- ifelse(grepl("TC", concordance_values$Model), "orange", "skyblue")

concordance_values$Label <- ifelse(grepl("TC", concordance_values$Model), "Includes TC", "Does Not Include TC")


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



tmp <- clinical_data_clean[clinical_data_clean$Sarculator_5_Year_OS == "Low", ]

tmp2 <- clinical_data_clean[clinical_data_clean$Sarculator_5_Year_OS == "High", ]


colnames(tmp)[colnames(tmp) == "Transcriptomic Clusters"] <- "Transcriptomic_Clusters"
tmp$Transcriptomic_Clusters <- ifelse(tmp$Transcriptomic_Clusters == "C1","C1", "Not C1")

colnames(tmp2)[colnames(tmp2) == "Transcriptomic Clusters"] <- "Transcriptomic_Clusters"
tmp2$Transcriptomic_Clusters <- ifelse(tmp2$Transcriptomic_Clusters == "C1","C1", "Not C1")



km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ Transcriptomic_Clusters, data=tmp)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival Sarculator Score <= 0.60")
x

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ Transcriptomic_Clusters, data=tmp2)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival Sarculator Score > 0.60")

x



