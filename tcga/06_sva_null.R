#!/usr/bin/env Rscript
# Do surrogate variables (SVA) remove the false positives that DESeq2 and edgeR make on the real-data permutation null?
#
# Null part: the SAME pseudo-mutant draws as tcga/04_permutation_null.R (same seeds, same 150 controls, same filter), so
# each draw can be compared with its own result without surrogate variables (results/tcga/perm/BRCA1_mut_iterNN.csv).
# For each draw: estimate the number of surrogate variables (SVs) with num.sv, estimate them with svaseq (protecting the
# case/control label and keeping subtype in both models), add them to the design, and rerun DESeq2, edgeR (QL),
# limma-voom and the rank-sum test. The SVs are correlated with marker genes (MKI67 proliferation, PTPRC immune,
# KRT19 epithelial) to give them biological meaning.
# Real part: the real BRCA1 mutants vs stratified controls, with SVs, for comparison with the null-with-SV range.
#
# Caveat: with real labels, SVs can absorb true signal that is correlated with them; with 11 mutants that is a real risk.
#
# Usage: Rscript tcga/06_sva_null.R [n_iter=10] [n_real=3] [n_ctrl=150]     (n_real = 0 skips the real part)
# Written without an R runtime available: run with `1 0` first. Resumable.

suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(edgeR); library(limma); library(sva)})
args   <- commandArgs(TRUE)
n_iter <- if (length(args) >= 1) as.integer(args[1]) else 10
n_real <- if (length(args) >= 2) as.integer(args[2]) else 3
n_ctrl <- if (length(args) >= 3) as.integer(args[3]) else 150
target <- "BRCA1_mut"
dir.create("results/tcga/sva", recursive = TRUE, showWarnings = FALSE)

se <- readRDS("data/tcga/brca_se.rds"); se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
cd <- as.data.frame(colData(se)); cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
counts_all <- assay(se, "unstranded"); gene_name <- rowData(se)$gene_name
case_n  <- table(cd$subtype[cd[[target]]])
pool    <- which(!cd$BRCA_mut)
markers <- c(proliferation = "MKI67", immune = "PTPRC", epithelial = "KRT19")
secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

strat_controls <- function(rest, n) {
  per <- split(rest, cd$subtype[rest]); per <- per[lengths(per) > 0]
  unname(unlist(lapply(per, function(ix) ix[sample.int(length(ix), min(length(ix), ceiling(n * length(ix) / length(rest))))])))
}
make_m <- function(idx, n_cases) data.frame(subtype = droplevels(cd$subtype[idx]),
  grp = factor(rep(c("case", "ctrl"), c(n_cases, length(idx) - n_cases)), levels = c("ctrl", "case")))

metric_names <- c("deseq2", "edger", "limmavoom", "wilcoxon", "p05_deseq2", "p05_edger", "p05_limmavoom", "p05_wilcoxon",
                  "sec_sva", "sec_deseq2", "sec_edger", "sec_limmavoom", "sec_wilcoxon", paste0("rho_", names(markers)))

analyze <- function(idx, m) {
  cm <- counts_all[, idx]
  keep <- filterByExpr(DGEList(cm), design = model.matrix(~ subtype + grp, m)); cm <- cm[keep, ]   # as in 04
  gn <- gene_name[keep]
  mod <- model.matrix(~ subtype + grp, m); mod0 <- model.matrix(~ subtype, m)
  out <- setNames(as.list(rep(NA_real_, length(metric_names))), metric_names)

  t0 <- Sys.time()
  n_sv <- min(10, num.sv(log(cm + 1), mod, vfilter = 10000, B = 20, seed = 1))
  base <- data.frame(n_case = sum(m$grp == "case"), n_ctrl = sum(m$grp == "ctrl"), genes_tested = nrow(cm), n_sv = n_sv)
  if (n_sv == 0) return(cbind(base, as.data.frame(out)))                       # nothing to adjust for
  sv <- svaseq(cm, mod, mod0, n.sv = n_sv, vfilter = 10000)$sv
  colnames(sv) <- paste0("SV", seq_len(n_sv)); m2 <- cbind(m, sv)
  out$sec_sva <- secs(t0)

  y  <- calcNormFactors(DGEList(cm)); lc <- cpm(y, log = TRUE, prior.count = 1)
  for (nm in names(markers)) {
    g <- which(gn == markers[[nm]])[1]
    if (!is.na(g)) out[[paste0("rho_", nm)]] <- max(abs(cor(sv, lc[g, ], method = "spearman")))   # strongest |Spearman rho| over SVs
  }

  fml <- as.formula(paste("~ subtype +", paste(c(colnames(sv), "grp"), collapse = " + ")))
  des <- model.matrix(fml, m2)
  calls <- function(q) sum(!is.na(q) & q < 0.05); p05 <- function(p) mean(p[!is.na(p)] < 0.05)

  t0 <- Sys.time()
  dds <- DESeq(DESeqDataSetFromMatrix(cm, m2, fml), quiet = TRUE)
  d <- results(dds, contrast = c("grp", "case", "ctrl")); rm(dds)
  out$sec_deseq2 <- secs(t0); out$deseq2 <- calls(d$padj); out$p05_deseq2 <- p05(d$pvalue)

  t0 <- Sys.time()
  yd <- estimateDisp(y, des); e <- glmQLFTest(glmQLFit(yd, des), coef = "grpcase")$table
  out$sec_edger <- secs(t0); out$edger <- calls(p.adjust(e$PValue, "BH")); out$p05_edger <- p05(e$PValue)

  t0 <- Sys.time()
  l <- topTable(eBayes(lmFit(voom(yd, des), des)), coef = "grpcase", number = Inf, sort.by = "none")
  out$sec_limmavoom <- secs(t0); out$limmavoom <- calls(p.adjust(l$P.Value, "BH")); out$p05_limmavoom <- p05(l$P.Value)

  t0 <- Sys.time()
  rs <- t(lm.fit(model.matrix(as.formula(paste("~ subtype +", paste(colnames(sv), collapse = " + "))), m2), t(lc))$residuals)
  is_case <- m$grp == "case"; n1 <- sum(is_case); n2 <- sum(!is_case)
  rk <- t(apply(rs, 1, rank)); U <- rowSums(rk[, is_case, drop = FALSE]) - n1 * (n1 + 1) / 2
  pw <- 2 * pnorm(-abs((U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
  out$sec_wilcoxon <- secs(t0); out$wilcoxon <- calls(p.adjust(pw, "BH")); out$p05_wilcoxon <- p05(pw)
  cbind(base, as.data.frame(out))
}

report <- function(kind, i, r) message(kind, " ", i, ": n_sv=", r$n_sv, " | calls deseq2/edger/voom/wilcoxon: ",
  paste(unlist(r[1, c("deseq2", "edger", "limmavoom", "wilcoxon")]), collapse = "/"),
  " | seconds sva/deseq2/edger: ", paste(round(unlist(r[1, c("sec_sva", "sec_deseq2", "sec_edger")])), collapse = "/"))

## Null draws (same construction and seeds as 04_permutation_null.R)
for (i in seq_len(n_iter)) {
  f <- sprintf("results/tcga/sva/%s_null_iter%02d.csv", target, i)
  if (file.exists(f)) { message("null ", i, " cached"); next }
  set.seed(1000 + i)
  cases <- unlist(lapply(names(case_n), function(s) {
    if (case_n[[s]] == 0) return(NULL)
    p <- pool[cd$subtype[pool] == s]; stopifnot(length(p) >= case_n[[s]])
    p[sample.int(length(p), case_n[[s]])]
  }))
  idx <- unname(c(cases, strat_controls(setdiff(pool, cases), n_ctrl)))
  r <- analyze(idx, make_m(idx, length(cases))); r$iter <- i
  write.csv(r, f, row.names = FALSE); report("null", i, r); gc()
}

## Real BRCA1 mutants vs stratified controls, with SVs
cases_real <- which(cd[[target]])
for (i in seq_len(n_real)) {
  f <- sprintf("results/tcga/sva/%s_real_iter%02d.csv", target, i)
  if (file.exists(f)) { message("real ", i, " cached"); next }
  set.seed(9000 + i)
  idx <- unname(c(cases_real, strat_controls(setdiff(pool, cases_real), n_ctrl)))
  r <- analyze(idx, make_m(idx, length(cases_real))); r$iter <- i
  write.csv(r, f, row.names = FALSE); report("real", i, r); gc()
}

## Summary: paired comparison with the null without SVs
read_all <- function(dir, pat) {
  f <- list.files(dir, pattern = pat, full.names = TRUE)
  x <- Filter(function(d) "wilcoxon" %in% names(d), lapply(f, read.csv))
  if (length(x)) do.call(rbind, x) else NULL
}
nul <- read_all("results/tcga/sva", paste0("^", target, "_null_iter[0-9]+\\.csv$"))
bas <- read_all("results/tcga/perm", paste0("^", target, "_iter[0-9]+\\.csv$"))
if (!is.null(bas)) bas <- bas[bas$design == "adjusted", ]
if (!is.null(nul) && !is.null(bas)) {
  rows <- lapply(c("deseq2", "edger", "limmavoom", "wilcoxon"), function(t) {
    j <- merge(bas[, c("iter", t, paste0("p05_", t))], nul[, c("iter", t, paste0("p05_", t))], by = "iter", suffixes = c("_base", "_sva"))
    data.frame(tool = t, n_paired = nrow(j),
               mean_calls_without_sva = mean(j[[paste0(t, "_base")]], na.rm = TRUE), mean_calls_with_sva = mean(j[[paste0(t, "_sva")]], na.rm = TRUE),
               mean_pct_p05_without_sva = 100 * mean(j[[paste0("p05_", t, "_base")]], na.rm = TRUE),
               mean_pct_p05_with_sva = 100 * mean(j[[paste0("p05_", t, "_sva")]], na.rm = TRUE))
  })
  summ <- do.call(rbind, rows)
  write.csv(summ, "results/tcga/sva_null_summary.csv", row.names = FALSE)
  cat("\nNull draws, paired with and without surrogate variables (subtype-adjusted in both):\n"); print(summ, row.names = FALSE)
  cat("\nSurrogate variables per draw:", paste(nul$n_sv, collapse = ", "), "\n")
  cat("Strongest |Spearman rho| of any SV with marker genes (mean over draws):\n"); print(round(colMeans(nul[, grep("^rho_", names(nul))], na.rm = TRUE), 2))
}
rl <- read_all("results/tcga/sva", paste0("^", target, "_real_iter[0-9]+\\.csv$"))
if (!is.null(rl)) { cat("\nReal BRCA1 mutants with SVs (one row per control draw):\n"); print(rl[, c("iter", "n_sv", "deseq2", "edger", "limmavoom", "wilcoxon")], row.names = FALSE) }
