#!/usr/bin/env Rscript
# Wilcoxon rank-sum check on TCGA-BRCA: does it have POWER (positive control) and is it CALIBRATED (null),
# and does it find anything in the REAL BRCA mutants beyond what the null produces?
#
# Test: rank-sum on log-CPM (normal approximation). "adjusted" regresses subtype out of log-CPM first and tests
# the residuals (an approximation to a subtype-stratified test). Fast: a few seconds per analysis.
#
# Sections:
#   1. Positive control: 100 Basal vs 100 LumA (should call MANY genes -> the test has power on a strong effect).
#   2. Null: pseudo-mutants drawn from wild-type tumors with the BRCA1 subtype mix (should call ~0 genes).
#   3. Real BRCA1 / BRCA2 / pooled mutants vs stratified controls, several control draws, naive and adjusted.
#
# Usage: Rscript tcga/05_wilcoxon_check.R [n_null=50] [n_draws=10] [n_ctrl=150]
# Written without an R runtime available: if it errors, paste the message. Uses ~1 GB; check `free -m` first
# if the permutation job is still running.

suppressPackageStartupMessages({library(SummarizedExperiment); library(edgeR)})
args   <- commandArgs(TRUE)
n_null  <- if (length(args) >= 1) as.integer(args[1]) else 50
n_draws <- if (length(args) >= 2) as.integer(args[2]) else 10
n_ctrl  <- if (length(args) >= 3) as.integer(args[3]) else 150
dir.create("results/tcga", recursive = TRUE, showWarnings = FALSE)

se <- readRDS("data/tcga/brca_se.rds")
se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
cd <- as.data.frame(colData(se)); cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
counts_all <- assay(se, "unstranded"); gene_name <- rowData(se)$gene_name
pool <- which(!cd$BRCA_mut)                                  # wild-type for BRCA1 and BRCA2

wilcox_rs <- function(idx, is_case, adjust) {
  cm <- counts_all[, idx]
  m <- data.frame(subtype = droplevels(cd$subtype[idx]), grp = factor(ifelse(is_case, "case", "ctrl"), levels = c("ctrl", "case")))
  keep <- filterByExpr(DGEList(cm), design = model.matrix(~ grp, m)); cm <- cm[keep, ]
  y <- calcNormFactors(DGEList(cm)); lc <- cpm(y, log = TRUE, prior.count = 1)       # genes x samples
  if (adjust && nlevels(m$subtype) > 1) lc <- t(lm.fit(model.matrix(~ subtype, m), t(lc))$residuals)
  n1 <- sum(is_case); n2 <- sum(!is_case)
  rk <- t(apply(lc, 1, rank))
  U  <- rowSums(rk[, is_case, drop = FALSE]) - n1 * (n1 + 1) / 2
  z  <- (U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  p  <- 2 * pnorm(-abs(z))
  data.frame(gene_id = rownames(cm), gene_name = gene_name[keep], p = p, q = p.adjust(p, "BH"))
}

strat_controls <- function(rest, n) {
  per <- split(rest, cd$subtype[rest]); per <- per[lengths(per) > 0]
  unname(unlist(lapply(per, function(ix) ix[sample.int(length(ix), min(length(ix), ceiling(n * length(ix) / length(rest))))])))
}
calls <- function(r) sum(r$q < 0.05, na.rm = TRUE)
summ_row <- function(analysis, design, v) data.frame(analysis = analysis, design = design, n_draws = length(v),
  mean_calls = mean(v), min_calls = min(v), max_calls = max(v), frac_draws_with_any_call = mean(v > 0))
out <- list()

## 1. Positive control
set.seed(1)
bas <- which(cd$subtype == "Basal"); lum <- which(cd$subtype == "LumA")
idx <- c(sample(bas, 100), sample(lum, 100))
r <- wilcox_rs(idx, cd$subtype[idx] == "Basal", adjust = FALSE)
cat("Positive control, Basal vs LumA (100 vs 100):", calls(r), "genes at q<0.05 of", nrow(r), "tested;",
    round(mean(r$p < 0.05) * 100, 1), "% raw p<0.05\n")
out[[1]] <- summ_row("positive_control_basal_vs_lumA", "naive", calls(r))

## 2. Null: pseudo-mutants with the real BRCA1 subtype mix
case_n <- table(cd$subtype[cd$BRCA1_mut]); v <- integer(n_null); p05 <- numeric(n_null)
for (i in seq_len(n_null)) {
  set.seed(5000 + i)
  cases <- unlist(lapply(names(case_n), function(s) if (case_n[[s]] == 0) NULL else {
    p <- pool[cd$subtype[pool] == s]; p[sample.int(length(p), case_n[[s]])] }))
  idx <- c(cases, strat_controls(setdiff(pool, cases), n_ctrl))
  r <- wilcox_rs(idx, seq_along(idx) <= length(cases), adjust = TRUE)
  v[i] <- calls(r); p05[i] <- mean(r$p < 0.05)
}
cat("Null (", n_null, " draws, adjusted): calls per draw mean ", mean(v), ", max ", max(v),
    "; mean raw p<0.05 = ", round(mean(p05) * 100, 1), "%\n", sep = "")
out[[2]] <- summ_row("null_pseudo_BRCA1", "adjusted", v)

## 3. Real mutants
for (target in c("BRCA1_mut", "BRCA2_mut", "BRCA_mut")) for (design in c("naive", "adjusted")) {
  v <- integer(n_draws); called <- list(); cases <- which(cd[[target]])
  for (i in seq_len(n_draws)) {
    set.seed(9000 + i)
    idx <- c(cases, strat_controls(setdiff(pool, cases), n_ctrl))
    r <- wilcox_rs(idx, seq_along(idx) <= length(cases), adjust = design == "adjusted")
    v[i] <- calls(r); called[[i]] <- r$gene_id[r$q < 0.05 & !is.na(r$q)]
  }
  cat(target, design, "- real mutants (n =", length(cases), "): calls per control draw:", paste(v, collapse = ", "), "\n")
  out[[length(out) + 1]] <- summ_row(target, design, v)
  tab <- sort(table(unlist(called)), decreasing = TRUE)
  if (length(tab)) {
    g <- data.frame(gene_id = names(tab), n_draws_called = as.integer(tab))
    g$gene_name <- gene_name[match(g$gene_id, rownames(se))]
    write.csv(g, sprintf("results/tcga/wilcoxon_called_%s_%s.csv", target, design), row.names = FALSE)
  }
}

res <- do.call(rbind, out); print(res, row.names = FALSE)
write.csv(res, "results/tcga/wilcoxon_check_summary.csv", row.names = FALSE)
