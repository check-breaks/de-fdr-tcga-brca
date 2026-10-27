#!/usr/bin/env Rscript
# Download open-access TCGA-BRCA data from the GDC and build analysis-ready objects.
# Written without network/R access; column names in GDC outputs change between releases -- check each step.
suppressPackageStartupMessages({library(TCGAbiolinks); library(SummarizedExperiment); library(dplyr)})
dir.create("data/tcga", recursive = TRUE, showWarnings = FALSE)

## 1. RNA-seq raw counts (STAR - Counts)
q <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts",
              sample.type = "Primary Tumor")
GDCdownload(q, directory = "data/tcga/gdc")
se <- GDCprepare(q, directory = "data/tcga/gdc")
# one tumor sample per patient
se <- se[, !duplicated(substr(colnames(se), 1, 12))]
colnames(se) <- substr(colnames(se), 1, 12)

## 2. Somatic mutations (open access MAF) -> BRCA1/BRCA2 status
qm <- GDCquery(project = "TCGA-BRCA", data.category = "Simple Nucleotide Variation",
               data.type = "Masked Somatic Mutation")
GDCdownload(qm, directory = "data/tcga/gdc")
maf <- GDCprepare(qm, directory = "data/tcga/gdc")
damaging <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site")
mut <- maf %>% filter(Hugo_Symbol %in% c("BRCA1", "BRCA2"), Variant_Classification %in% damaging) %>%
  mutate(patient = substr(Tumor_Sample_Barcode, 1, 12)) %>% distinct(patient, Hugo_Symbol)

cd <- as.data.frame(colData(se))
cd$patient <- rownames(cd)
cd$BRCA1_mut <- cd$patient %in% mut$patient[mut$Hugo_Symbol == "BRCA1"]
cd$BRCA2_mut <- cd$patient %in% mut$patient[mut$Hugo_Symbol == "BRCA2"]
cd$BRCA_mut  <- cd$BRCA1_mut | cd$BRCA2_mut

## 3. PAM50 subtype (the key confounder: BRCA1-mutant tumors are mostly Basal-like)
sub <- TCGAquery_subtype("BRCA")
pam <- sub[, c("patient", "BRCA_Subtype_PAM50")]    # verify column name in your TCGAbiolinks version
cd <- left_join(cd, pam, by = "patient"); rownames(cd) <- cd$patient
colData(se) <- DataFrame(cd[colnames(se), ])

saveRDS(se, "data/tcga/brca_se.rds")
print(table(PAM50 = cd$BRCA_Subtype_PAM50, BRCA1 = cd$BRCA1_mut, useNA = "ifany"))
cat("BRCA1 mut:", sum(cd$BRCA1_mut), " BRCA2 mut:", sum(cd$BRCA2_mut), "\n")
