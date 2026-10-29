#!/usr/bin/env Rscript
# Improved Kaplan-Meier plot for BRCA1 expression: time in years, x-axis truncated at 12 years (few patients
# remain at risk beyond that, so late steps are driven by a handful of people), risk table included.
# Written without an R runtime available. Usage: Rscript tcga/03b_km_plot.R
suppressPackageStartupMessages({library(SummarizedExperiment); library(survival); library(survminer)})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
se <- readRDS("data/tcga/brca_se.rds"); cd <- as.data.frame(colData(se))
lcpm <- log2(edgeR::cpm(assay(se, "unstranded")) + 1)
cd$BRCA1_expr <- lcpm[which(rowData(se)$gene_name == "BRCA1")[1], ]
cd$time  <- ifelse(!is.na(cd$days_to_death), cd$days_to_death, cd$days_to_last_follow_up) / 365.25
cd$event <- as.integer(cd$vital_status == "Dead")
cd <- cd[!is.na(cd$time) & !is.na(cd$event), ]
cd$BRCA1 <- factor(ifelse(cd$BRCA1_expr > median(cd$BRCA1_expr, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
fit <- survfit(Surv(time, event) ~ BRCA1, data = cd)
p <- ggsurvplot(fit, data = cd, pval = TRUE, pval.method = TRUE, risk.table = TRUE, xlim = c(0, 12), break.time.by = 2,
                xlab = "Time (years)", ylab = "Overall survival", legend.title = "BRCA1",
                legend.labs = c("Low", "High"), risk.table.height = 0.25,
                title = sprintf("BRCA1 expression (median split): n = %d, %d deaths", sum(fit$n), sum(fit$n.event)))
message("patients excluded from the KM fit: ", nrow(cd) - sum(fit$n))
png("results/figures/km_brca1_v2.png", width = 1700, height = 1200, res = 200); print(p); dev.off()
cat("saved results/figures/km_brca1_v2.png\n")
