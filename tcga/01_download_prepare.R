#!/usr/bin/env Rscript
# Download open-access TCGA-BRCA data from the GDC and build the analysis object.
#   - Memory-lean: reads only the raw-count column of each STAR file (GDCprepare on all assays needs several GB).
#   - Clinical data come from the GDC API. Follow-up days are stored in `follow_ups` (not `diagnoses`), so both are queried;
#     without this, living patients have no survival time and are silently dropped from survival models.
# Usage: Rscript tcga/01_download_prepare.R [output.rds]      (default output: data/tcga/brca_se.rds)
# GDC field names change between releases: check the printouts at each step.
options(timeout = 100000)
suppressPackageStartupMessages({library(TCGAbiolinks); library(SummarizedExperiment)
  library(dplyr); library(data.table)})
args <- commandArgs(TRUE)
out_file <- if (length(args)) args[1] else "data/tcga/brca_se.rds"
dir.create("data/tcga", recursive = TRUE, showWarnings = FALSE)
gdc <- "data/tcga/gdc"

## 1. Raw counts (one primary tumor per patient)
q <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts",
              sample.type = "Primary Tumor")
GDCdownload(q, directory = gdc, files.per.chunk = 10)
res <- getResults(q)
stopifnot(all(grepl("^TCGA-", res$cases)))
files <- list.files(file.path(gdc, "TCGA-BRCA"), pattern = "star_gene_counts.tsv$",
                    recursive = TRUE, full.names = TRUE)
res$path <- setNames(files, basename(files))[res$file_name]
stopifnot(!anyNA(res$path))
res$patient <- substr(res$cases, 1, 12)
res <- res[!duplicated(res$patient), ]
message("samples: ", nrow(res))

first <- fread(res$path[1], skip = "gene_id", select = c("gene_id", "gene_name", "gene_type"))
first <- first[!grepl("^N_", gene_id)]
mat <- matrix(0L, nrow(first), nrow(res), dimnames = list(first$gene_id, res$patient))
for (i in seq_len(nrow(res))) {
  d <- fread(res$path[i], skip = "gene_id", select = c("gene_id", "unstranded"))
  d <- d[!grepl("^N_", gene_id)]
  stopifnot(identical(d$gene_id, first$gene_id))
  mat[, i] <- as.integer(d$unstranded)
  if (i %% 100 == 0) message(i, " / ", nrow(res))
}

## 2. Clinical via the GDC API (GDCquery_clinic has a data.table bug in some TCGAbiolinks versions)
filt <- list(op = "in", content = list(field = "project.project_id", value = list("TCGA-BRCA")))
fields <- c("submitter_id", "demographic.vital_status", "demographic.days_to_death",
            "diagnoses.days_to_last_follow_up", "diagnoses.age_at_diagnosis",
            "diagnoses.ajcc_pathologic_stage", "follow_ups.days_to_follow_up")
r <- httr::GET("https://api.gdc.cancer.gov/cases",
               query = list(filters = jsonlite::toJSON(filt, auto_unbox = TRUE),
                            fields = paste(fields, collapse = ","), size = "2000", format = "JSON"))
httr::stop_for_status(r)
h <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), flatten = TRUE)$data$hits
message("clinical columns returned: ", paste(names(h), collapse = ", "))

# first diagnosis row per case (age, stage, and a follow-up value where present)
dx <- lapply(h$diagnoses, function(x) if (is.data.frame(x) && nrow(x)) x[1, , drop = FALSE] else NULL)
gx <- function(n) vapply(dx, function(x) if (is.null(x) || !n %in% names(x) || is.null(x[[n]])) NA_character_ else as.character(x[[n]]), character(1))
gh <- function(n) if (n %in% names(h)) h[[n]] else rep(NA, nrow(h))
# latest value over all rows of a list-column of data frames (a case can have several follow-ups)
mx <- function(name, col) {
  if (!name %in% names(h)) return(rep(NA_real_, nrow(h)))
  vapply(h[[name]], function(x) {
    if (is.data.frame(x) && col %in% names(x)) {
      v <- suppressWarnings(as.numeric(x[[col]]))
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    } else NA_real_
  }, numeric(1))
}
dx_fu <- as.numeric(gx("days_to_last_follow_up"))
last_fu <- ifelse(is.na(dx_fu),
                  pmax(mx("diagnoses", "days_to_last_follow_up"), mx("follow_ups", "days_to_follow_up"), na.rm = TRUE),
                  dx_fu)
clin <- data.frame(submitter_id = h$submitter_id,
                   vital_status = gh("demographic.vital_status"),
                   days_to_death = as.numeric(gh("demographic.days_to_death")),
                   days_to_last_follow_up = last_fu,
                   age_at_diagnosis = as.numeric(gx("age_at_diagnosis")),
                   ajcc_pathologic_stage = gx("ajcc_pathologic_stage"),
                   stringsAsFactors = FALSE)
clin <- clin[!duplicated(clin$submitter_id), ]
message("clinical rows: ", nrow(clin))

## 3. Somatic mutations -> BRCA1/BRCA2 status
qm <- GDCquery(project = "TCGA-BRCA", data.category = "Simple Nucleotide Variation",
               data.type = "Masked Somatic Mutation")
GDCdownload(qm, directory = gdc, files.per.chunk = 10)
maf <- GDCprepare(qm, directory = gdc)
damaging <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site")
mut <- maf %>% filter(Hugo_Symbol %in% c("BRCA1", "BRCA2"), Variant_Classification %in% damaging) %>%
  mutate(patient = substr(Tumor_Sample_Barcode, 1, 12)) %>% distinct(patient, Hugo_Symbol)
rm(maf); gc()

## 4. PAM50 subtype (key confounder)
sub <- TCGAquery_subtype("BRCA")
if (!"BRCA_Subtype_PAM50" %in% colnames(sub)) stop("PAM50 column not found; columns are: ", paste(colnames(sub), collapse = ", "))
pam <- sub[, c("patient", "BRCA_Subtype_PAM50")]

cd <- data.frame(patient = colnames(mat)) %>%
  left_join(clin, by = c("patient" = "submitter_id")) %>%
  left_join(pam, by = "patient") %>%
  mutate(BRCA1_mut = patient %in% mut$patient[mut$Hugo_Symbol == "BRCA1"],
         BRCA2_mut = patient %in% mut$patient[mut$Hugo_Symbol == "BRCA2"],
         BRCA_mut = BRCA1_mut | BRCA2_mut)

## Guard against silent loss of survival data (the earlier version dropped every living patient)
no_time <- is.na(cd$days_to_death) & is.na(cd$days_to_last_follow_up)
cat("patients with no usable survival time:", sum(no_time), "of", nrow(cd), "\n")
print(table(vital_status = cd$vital_status, no_time = no_time, useNA = "ifany"))
if (mean(no_time) > 0.05) warning("More than 5% of patients have no usable survival time; check the GDC follow-up fields")

se <- SummarizedExperiment(assays = list(unstranded = mat),
                           rowData = DataFrame(as.data.frame(first)),
                           colData = DataFrame(cd, row.names = cd$patient))
saveRDS(se, out_file)
cat("saved", out_file, "\n")
print(table(PAM50 = cd$BRCA_Subtype_PAM50, BRCA1 = cd$BRCA1_mut, useNA = "ifany"))
cat("BRCA1 mut:", sum(cd$BRCA1_mut), " BRCA2 mut:", sum(cd$BRCA2_mut), "\n")
