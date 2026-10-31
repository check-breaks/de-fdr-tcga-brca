#!/usr/bin/env Rscript
# Permutation null within ONE PAM50 subtype (wild-type tumors only; pseudo-cases and controls both from that subtype).
# Question: does the DESeq2/edgeR inflation disappear when subtype is constant?
# Usage: Rscript tcga/04d_permutation_single_subtype.R [n_iter=3] [subtype=LumA] [n_case=11] [n_ctrl=150]
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(edgeR); library(limma)})
args    <- commandArgs(TRUE)
n_iter  <- if (length(args) >= 1) as.integer(args[1]) else 3
sub_t   <- if (length(args) >= 2) args[2] else "LumA"
n_case  <- if (length(args) >= 3) as.integer(args[3]) else 11
n_ctrl  <- if (length(args) >= 4) as.integer(args[4]) else 150
outdir  <- "results/tcga/perm_single"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

se <- readRDS("data/tcga/brca_se.rds")
se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
cd <- as.data.frame(colData(se))
counts_all <- assay(se, "unstranded")
pool <- which(cd$BRCA_Subtype_PAM50 == sub_t & cd$BRCA_mut %in% FALSE)
cat("Subtype:", sub_t, "| wild-type pool:", length(pool), "| cases:", n_case, "| controls:", n_ctrl, "\n")
stopifnot(length(pool) >= n_case + n_ctrl)

secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))
calls <- function(q) sum(!is.na(q) & q < 0.05)
p05 <- function(p) mean(p[!is.na(p)] < 0.05)

for (i in seq_len(n_iter)) {
  f <- sprintf("%s/%s_c%d_n%d_iter%02d.csv", outdir, sub_t, n_case, n_ctrl, i)
  if (file.exists(f)) { message("iter ", i, " cached"); next }
  set.seed(2000 + i)
  idx <- sample(pool, n_case + n_ctrl)
  m <- data.frame(grp = factor(rep(c("case", "ctrl"), c(n_case, n_ctrl)), levels = c("ctrl", "case")))
  des <- model.matrix(~ grp, m)
  cm <- counts_all[, idx]
  cm <- cm[filterByExpr(DGEList(cm), design = des), ]

  t0 <- Sys.time()
  dds <- DESeq(DESeqDataSetFromMatrix(cm, m, ~ grp), quiet = TRUE)
  d <- results(dds, contrast = c("grp", "case", "ctrl")); s1 <- secs(t0); rm(dds)
  t0 <- Sys.time()
  y <- calcNormFactors(DGEList(cm)); yd <- estimateDisp(y, des)
  e <- glmQLFTest(glmQLFit(yd, des), coef = "grpcase")$table; s2 <- secs(t0)
  t0 <- Sys.time()
  l <- topTable(eBayes(lmFit(voom(yd, des), des)), coef = "grpcase", number = Inf, sort.by = "none"); s3 <- secs(t0)
  lc <- cpm(y, log = TRUE, prior.count = 1)
  is_case <- m$grp == "case"; n1 <- sum(is_case); n2 <- sum(!is_case)
  rk <- t(apply(lc, 1, rank))
  U <- rowSums(rk[, is_case, drop = FALSE]) - n1 * (n1 + 1) / 2
  pw <- 2 * pnorm(-abs((U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))

  out <- data.frame(subtype = sub_t, n_case = n_case, n_ctrl = n_ctrl, genes_tested = nrow(cm),
                    deseq2 = calls(d$padj), edger = calls(p.adjust(e$PValue, "BH")),
                    limmavoom = calls(p.adjust(l$P.Value, "BH")), wilcoxon = calls(p.adjust(pw, "BH")),
                    p05_deseq2 = p05(d$pvalue), p05_edger = p05(e$PValue),
                    p05_limmavoom = p05(l$P.Value), p05_wilcoxon = p05(pw),
                    sec_deseq2 = s1, sec_edger = s2, iter = i)
  write.csv(out, f, row.names = FALSE)
  message("iter ", i, " calls deseq2/edger/voom/wilcoxon: ", paste(out[1, 5:8], collapse = "/"),
          " | raw p<0.05 deseq2/edger/voom/wilcoxon: ", paste(round(unlist(out[1, 9:12]), 3), collapse = "/"))
  gc()
}

files <- list.files(outdir, pattern = "iter[0-9]+\\.csv$", full.names = TRUE)
all <- do.call(rbind, lapply(files, read.csv))
agg <- aggregate(cbind(deseq2, edger, limmavoom, wilcoxon, p05_deseq2, p05_edger, p05_limmavoom, p05_wilcoxon) ~
                   subtype + n_case + n_ctrl, data = all, FUN = mean)
agg$n_iter <- aggregate(iter ~ subtype + n_case + n_ctrl, data = all, FUN = length)$iter
write.csv(agg, "results/tcga/single_subtype_null_summary.csv", row.names = FALSE)
print(agg, row.names = FALSE)
