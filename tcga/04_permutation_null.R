#!/usr/bin/env Rscript
# Permutation null on REAL TCGA-BRCA data (v2: adds Wilcoxon, per-tool timing, n_ctrl argument).
#
# Replace the real mutants with "pseudo-mutants" drawn at random from wild-type tumors, keeping the same subtype
# composition as the real mutants (BRCA1: 7 Basal / 2 LumA / 2 LumB). Calls on pseudo-mutants are almost all false
# positives, so this measures each method's false-positive behaviour in this unbalanced, heterogeneous design.
#
# Methods: DESeq2, edgeR (QL), limma-voom, and a Wilcoxon rank-sum test on log-CPM. For the "adjusted" design the
# Wilcoxon is run on log-CPM residuals after regressing out subtype (an approximation to a stratified test).
#
# Caveats: some "wild-type" tumors may be functionally BRCA-deficient; few iterations give rough rates.
#
# Usage: Rscript tcga/04_permutation_null.R [n_iter=10] [designs=adjusted|both] [target=BRCA1_mut] [n_ctrl=150]
# Written without an R runtime available. Resumable: finished iterations are skipped.

suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(edgeR); library(limma)})
args    <- commandArgs(TRUE)
n_iter  <- if (length(args) >= 1) as.integer(args[1]) else 10
designs <- if (length(args) >= 2 && args[2] == "both") c("adjusted", "naive") else "adjusted"
target  <- if (length(args) >= 3) args[3] else "BRCA1_mut"       # BRCA1_mut | BRCA2_mut | BRCA_mut
n_ctrl  <- if (length(args) >= 4) as.integer(args[4]) else 150
dir.create("results/tcga/perm", recursive = TRUE, showWarnings = FALSE)

se <- readRDS("data/tcga/brca_se.rds")
se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
cd <- as.data.frame(colData(se)); cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
counts_all <- assay(se, "unstranded")
stopifnot(target %in% colnames(cd))

case_n <- table(cd$subtype[cd[[target]]])          # subtype composition of the REAL mutants
pool   <- which(!cd$BRCA_mut)                      # wild-type for both BRCA1 and BRCA2
cat("Target:", target, "| controls:", n_ctrl, "| real mutant composition:\n"); print(case_n)
cat("Wild-type pool size:", length(pool), "| designs:", paste(designs, collapse = ", "), "\n")

secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

fit_one <- function(idx, m, design) {
  adj <- design == "adjusted" && nlevels(m$subtype) > 1
  fml <- if (adj) ~ subtype + grp else ~ grp
  des <- model.matrix(fml, m)
  cm  <- counts_all[, idx]
  keep <- filterByExpr(DGEList(cm), design = des); cm <- cm[keep, ]

  t0 <- Sys.time()
  dds <- DESeq(DESeqDataSetFromMatrix(cm, m, fml), quiet = TRUE)
  d <- results(dds, contrast = c("grp", "case", "ctrl"))            # padj includes independent filtering, as in 02
  s_deseq2 <- secs(t0); rm(dds)

  t0 <- Sys.time()
  y <- calcNormFactors(DGEList(cm)); yd <- estimateDisp(y, des)
  e <- glmQLFTest(glmQLFit(yd, des), coef = "grpcase")$table
  s_edger <- secs(t0)

  t0 <- Sys.time()
  l <- topTable(eBayes(lmFit(voom(yd, des), des)), coef = "grpcase", number = Inf, sort.by = "none")
  s_voom <- secs(t0)

  t0 <- Sys.time()
  lc <- cpm(y, log = TRUE, prior.count = 1)                          # genes x samples
  if (adj) lc <- t(lm.fit(model.matrix(~ subtype, m), t(lc))$residuals)
  is_case <- m$grp == "case"; n1 <- sum(is_case); n2 <- sum(!is_case)
  rk <- t(apply(lc, 1, rank))                                        # rank-sum test, normal approximation
  U  <- rowSums(rk[, is_case, drop = FALSE]) - n1 * (n1 + 1) / 2
  z  <- (U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  pw <- 2 * pnorm(-abs(z))
  s_wilcox <- secs(t0)

  calls <- function(q) sum(!is.na(q) & q < 0.05)
  p05   <- function(p) mean(p[!is.na(p)] < 0.05)                    # ~0.05 if calibrated under the null
  data.frame(design = design, n_case = sum(m$grp == "case"), n_ctrl = sum(m$grp == "ctrl"), genes_tested = nrow(cm),
             deseq2 = calls(d$padj), edger = calls(p.adjust(e$PValue, "BH")),
             limmavoom = calls(p.adjust(l$P.Value, "BH")), wilcoxon = calls(p.adjust(pw, "BH")),
             p05_deseq2 = p05(d$pvalue), p05_edger = p05(e$PValue), p05_limmavoom = p05(l$P.Value), p05_wilcoxon = p05(pw),
             sec_deseq2 = s_deseq2, sec_edger = s_edger, sec_limmavoom = s_voom, sec_wilcoxon = s_wilcox)
}

for (i in seq_len(n_iter)) {
  f <- sprintf("results/tcga/perm/%s_iter%02d.csv", target, i)
  if (file.exists(f)) { message("iter ", i, " cached"); next }
  set.seed(1000 + i)
  cases <- unlist(lapply(names(case_n), function(s) {
    if (case_n[[s]] == 0) return(NULL)
    p <- pool[cd$subtype[pool] == s]; stopifnot(length(p) >= case_n[[s]])
    p[sample.int(length(p), case_n[[s]])]
  }))
  rest <- setdiff(pool, cases)
  per  <- split(rest, cd$subtype[rest]); per <- per[lengths(per) > 0]
  ctrl <- unlist(lapply(per, function(ix)
    ix[sample.int(length(ix), min(length(ix), ceiling(n_ctrl * length(ix) / length(rest))))]))
  idx <- unname(c(cases, ctrl))
  m <- data.frame(subtype = droplevels(cd$subtype[idx]),
                  grp = factor(rep(c("case", "ctrl"), c(length(cases), length(ctrl))), levels = c("ctrl", "case")))
  out <- do.call(rbind, lapply(designs, function(dn) fit_one(idx, m, dn)))
  out$iter <- i
  write.csv(out, f, row.names = FALSE)
  message("iter ", i, " [", out$design[1], "] calls deseq2/edger/voom/wilcoxon: ",
          paste(out[1, c("deseq2", "edger", "limmavoom", "wilcoxon")], collapse = "/"),
          " | seconds: ", paste(round(unlist(out[1, c("sec_deseq2", "sec_edger", "sec_limmavoom", "sec_wilcoxon")])), collapse = "/"))
  gc()
}

## Summary across finished iterations (skips files from the older version that lack the Wilcoxon columns)
files <- list.files("results/tcga/perm", pattern = paste0("^", target, "_iter[0-9]+\\.csv$"), full.names = TRUE)
all <- Filter(function(x) "wilcoxon" %in% names(x), lapply(files, read.csv))
all <- do.call(rbind, all)
rows <- list()
for (dn in unique(all$design)) for (tool in c("deseq2", "edger", "limmavoom", "wilcoxon")) {
  a <- all[all$design == dn, ]
  rows[[length(rows) + 1]] <- data.frame(
    target = target, design = dn, tool = tool, n_iter = nrow(a),
    mean_calls = mean(a[[tool]]), median_calls = median(a[[tool]]), max_calls = max(a[[tool]]),
    frac_iters_with_any_call = mean(a[[tool]] > 0),
    mean_frac_raw_p_below_0.05 = mean(a[[paste0("p05_", tool)]]),
    mean_seconds = mean(a[[paste0("sec_", tool)]]))
}
summ <- do.call(rbind, rows)
write.csv(summ, sprintf("results/tcga/permutation_null_summary_%s.csv", target), row.names = FALSE)
print(summ, row.names = FALSE)
