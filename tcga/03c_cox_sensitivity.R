#!/usr/bin/env Rscript
# Sensitivity checks for the BRCA1 expression Cox model (see methods_background.md, section 9.5):
#   (a) all patients with a recorded time, (b) excluding the 21 patients with time <= 0,
#   (c) (b) with subtype as a stratification variable (each subtype gets its own baseline hazard).
# Written without an R runtime available. Usage: Rscript tcga/03c_cox_sensitivity.R
suppressPackageStartupMessages({library(SummarizedExperiment); library(survival)})
se <- readRDS("data/tcga/brca_se.rds"); cd <- as.data.frame(colData(se))
lcpm <- log2(edgeR::cpm(assay(se, "unstranded")) + 1)
cd$BRCA1_expr <- lcpm[which(rowData(se)$gene_name == "BRCA1")[1], ]
cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
cd$time  <- ifelse(!is.na(cd$days_to_death), cd$days_to_death, cd$days_to_last_follow_up) / 365.25
cd$event <- as.integer(cd$vital_status == "Dead")
cd$age   <- cd$age_at_diagnosis / 365.25
cat("patients with time <= 0:", sum(cd$time <= 0, na.rm = TRUE), "\n")

report <- function(label, fml, d) {
  r <- tryCatch({
    m <- coxph(fml, data = d); ci <- summary(m)$conf.int["BRCA1_expr", ]; p <- summary(m)$coefficients["BRCA1_expr", "Pr(>|z|)"]
    sprintf("n=%d events=%d  BRCA1 HR %.2f (%.2f to %.2f), p = %.3f", m$n, m$nevent, ci[1], ci[3], ci[4], p)
  }, error = function(e) paste("failed:", conditionMessage(e)))
  cat(sprintf("%-46s %s\n", label, r))
}
pos <- cd[!is.na(cd$time) & cd$time > 0, ]
report("(a) all patients with a time",            Surv(time, event) ~ BRCA1_expr + subtype + age, cd[!is.na(cd$time), ])
report("(b) excluding time <= 0",                 Surv(time, event) ~ BRCA1_expr + subtype + age, pos)
report("(c) excluding time <= 0, strata(subtype)", Surv(time, event) ~ BRCA1_expr + strata(subtype) + age, pos)
