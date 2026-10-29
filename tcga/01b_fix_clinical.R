#!/usr/bin/env Rscript
# Fix missing follow-up times for living patients in data/tcga/brca_se.rds.
# Queries the GDC API for follow-up days from BOTH diagnoses and follow_ups and takes the latest available.
# Backs up the object first. Written without network access: check the before/after printout.
suppressPackageStartupMessages({library(SummarizedExperiment); library(httr); library(jsonlite)})
options(timeout = 600)
file.copy("data/tcga/brca_se.rds", "data/tcga/brca_se_backup.rds", overwrite = FALSE)
se <- readRDS("data/tcga/brca_se.rds")

filt <- list(op = "in", content = list(field = "project.project_id", value = list("TCGA-BRCA")))
fields <- c("submitter_id", "diagnoses.days_to_last_follow_up", "follow_ups.days_to_follow_up")
r <- GET("https://api.gdc.cancer.gov/cases",
         query = list(filters = toJSON(filt, auto_unbox = TRUE), fields = paste(fields, collapse = ","),
                      size = "2000", format = "JSON"))
stop_for_status(r)
h <- fromJSON(content(r, as = "text", encoding = "UTF-8"), flatten = TRUE)$data$hits
cat("columns returned:", paste(names(h), collapse = ", "), "\n")

# latest value per case from a list-column of data frames
mx <- function(lst, col) vapply(lst, function(x) {
  if (is.data.frame(x) && col %in% names(x)) {
    v <- suppressWarnings(as.numeric(x[[col]]))
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  } else NA_real_
}, numeric(1))
dx <- if ("diagnoses" %in% names(h)) mx(h$diagnoses, "days_to_last_follow_up") else rep(NA_real_, nrow(h))
fu <- if ("follow_ups" %in% names(h)) mx(h$follow_ups, "days_to_follow_up") else rep(NA_real_, nrow(h))
cat("cases with a follow-up time: from diagnoses =", sum(!is.na(dx)), "| from follow_ups =", sum(!is.na(fu)), "\n")
tab <- data.frame(patient = h$submitter_id, last_fu = pmax(dx, fu, na.rm = TRUE))

cd <- as.data.frame(colData(se)); cd$patient <- rownames(cd)
cat("before: missing days_to_last_follow_up =", sum(is.na(cd$days_to_last_follow_up)), "of", nrow(cd), "\n")
new <- tab$last_fu[match(cd$patient, tab$patient)]
cd$days_to_last_follow_up <- ifelse(is.na(cd$days_to_last_follow_up), new, cd$days_to_last_follow_up)
cat("after:  missing days_to_last_follow_up =", sum(is.na(cd$days_to_last_follow_up)), "of", nrow(cd), "\n")
cat("patients with no usable time (neither death nor follow-up day), by vital status:\n")
print(table(vital = cd$vital_status, no_time = is.na(cd$days_to_last_follow_up) & is.na(cd$days_to_death), useNA = "ifany"))

colData(se) <- DataFrame(cd, row.names = rownames(cd))
saveRDS(se, "data/tcga/brca_se.rds")
cat("saved. Now rerun tcga/03_survival_expression.R\n")
