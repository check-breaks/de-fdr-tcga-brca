#!/usr/bin/env Rscript
# BRCA1/BRCA2 expression by subtype and survival (KM + Cox adjusted for subtype/age).
suppressPackageStartupMessages({library(SummarizedExperiment); library(survival); library(survminer); library(dplyr); library(ggplot2)})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
se <- readRDS("data/tcga/brca_se.rds"); cd <- as.data.frame(colData(se))
lcpm <- log2(edgeR::cpm(assay(se, "unstranded")) + 1)
sym <- rowData(se)$gene_name
for (g in c("BRCA1", "BRCA2")) cd[[paste0(g, "_expr")]] <- lcpm[which(sym == g)[1], ]
cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
ggplot(filter(cd, !is.na(subtype)), aes(subtype, BRCA1_expr, fill = subtype)) + geom_boxplot() + theme_bw() +
  labs(y = "BRCA1 log2 CPM", title = "BRCA1 expression by PAM50 subtype")
ggsave("results/figures/brca1_by_subtype.png", width = 6, height = 4, dpi = 150)

cd$time <- ifelse(!is.na(cd$days_to_death), cd$days_to_death, cd$days_to_last_follow_up) / 365.25
cd$event <- as.integer(cd$vital_status == "Dead")
cd$age <- cd$age_at_diagnosis / 365.25                       # verify units for your release
cd$BRCA1_high <- cd$BRCA1_expr > median(cd$BRCA1_expr, na.rm = TRUE)
fit <- survfit(Surv(time, event) ~ BRCA1_high, data = cd)
ggsave("results/figures/km_brca1.png", ggsurvplot(fit, data = cd, pval = TRUE, risk.table = TRUE)$plot, width = 6, height = 4.5, dpi = 150)
print(summary(coxph(Surv(time, event) ~ BRCA1_expr + subtype + age, data = cd)))
