#!/usr/bin/env Rscript
# BRCA-mutant vs wild-type DE on TCGA-BRCA: naive vs subtype-adjusted, DESeq2 / edgeR / limma-voom.
# Controls are subsampled (stratified by subtype) to keep memory low on an 8 GB machine.
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(edgeR); library(limma)})
dir.create("results/tcga", recursive = TRUE, showWarnings = FALSE)
se <- readRDS("data/tcga/brca_se.rds")
se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
cd <- as.data.frame(colData(se))
cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
counts_all <- assay(se, "unstranded")
gene_name <- rowData(se)$gene_name

cat("BRCA1 mutants by subtype:\n"); print(table(cd$subtype, cd$BRCA1_mut))
cat("BRCA2 mutants by subtype:\n"); print(table(cd$subtype, cd$BRCA2_mut))
cat("BRCA1 / BRCA2 overlap:\n"); print(table(BRCA1 = cd$BRCA1_mut, BRCA2 = cd$BRCA2_mut))

sample_controls <- function(is_case, n_ctrl = 300) {
  ctrl <- which(!is_case)
  per <- split(ctrl, cd$subtype[ctrl]); per <- per[lengths(per) > 0]
  keep <- unlist(lapply(per, function(ix)
    ix[sample.int(length(ix), min(length(ix), ceiling(n_ctrl * length(ix) / length(ctrl))))]))
  c(which(is_case), keep)
}

run_all <- function(idx, is_case, adjust, tag) {
  m <- cd[idx, c("subtype", "BRCA_Subtype_PAM50")]
  m$grp <- factor(ifelse(is_case[idx], "case", "ctrl"), levels = c("ctrl", "case"))
  m$subtype <- droplevels(m$subtype); m <- m[, c("subtype", "grp")]
  fml <- if (adjust && nlevels(m$subtype) > 1) ~ subtype + grp else ~ grp
  des <- model.matrix(fml, m)
  f <- sprintf("results/tcga/de_%s.csv", tag)
  if (file.exists(f)) {
    tab <- read.csv(f); sig <- function(q) !is.na(q) & q < 0.05
    s <- c(n_case = sum(m$grp == "case"), n_ctrl = sum(m$grp == "ctrl"), genes_tested = nrow(tab),
           deseq2 = sum(sig(tab$q_deseq2)), edger = sum(sig(tab$q_edger)),
           limmavoom = sum(sig(tab$q_limmavoom)),
           all_three = sum(sig(tab$q_deseq2) & sig(tab$q_edger) & sig(tab$q_limmavoom)))
    message(tag, " (cached): ", paste(names(s), s, collapse = "  "))
    return(s)
  }
  cm <- counts_all[, idx]
  keep <- filterByExpr(DGEList(cm), design = des)
  cm <- cm[keep, ]
  dds <- DESeq(DESeqDataSetFromMatrix(cm, m, fml), quiet = TRUE)
  d <- results(dds, contrast = c("grp", "case", "ctrl"))
  y <- calcNormFactors(DGEList(cm)); y <- estimateDisp(y, des)
  e <- glmQLFTest(glmQLFit(y, des), coef = "grpcase")$table
  l <- topTable(eBayes(lmFit(voom(y, des), des)), coef = "grpcase", number = Inf, sort.by = "none")
  tab <- data.frame(gene_id = rownames(cm), gene_name = gene_name[keep], lfc = d$log2FoldChange,
                    q_deseq2 = d$padj, q_edger = p.adjust(e$PValue, "BH"),
                    q_limmavoom = p.adjust(l$P.Value, "BH"))
  write.csv(tab, sprintf("results/tcga/de_%s.csv", tag), row.names = FALSE)
  sig <- function(q) !is.na(q) & q < 0.05
  s <- c(n_case = sum(m$grp == "case"), n_ctrl = sum(m$grp == "ctrl"), genes_tested = nrow(tab),
         deseq2 = sum(sig(tab$q_deseq2)), edger = sum(sig(tab$q_edger)),
         limmavoom = sum(sig(tab$q_limmavoom)),
         all_three = sum(sig(tab$q_deseq2) & sig(tab$q_edger) & sig(tab$q_limmavoom)))
  message(tag, ": ", paste(names(s), s, collapse = "  "))
  rm(dds); gc(); s
}

res <- list(); set.seed(1)
bas <- which(cd$subtype == "Basal"); lum <- which(cd$subtype == "LumA")
res$posctrl_basal_vs_lumA <- run_all(c(sample(bas, 100), sample(lum, 100)),
                                     cd$subtype == "Basal", FALSE, "posctrl_basal_vs_lumA")
for (g in c("BRCA1_mut", "BRCA2_mut", "BRCA_mut")) {
  is_case <- cd[[g]]; idx <- sample_controls(is_case)
  res[[paste0(g, "_naive")]]    <- run_all(idx, is_case, FALSE, paste0(g, "_naive"))
  res[[paste0(g, "_adjusted")]] <- run_all(idx, is_case, TRUE,  paste0(g, "_adjusted"))
}
res$BRCA1_basal_only <- run_all(bas, cd$BRCA1_mut, FALSE, "BRCA1_basal_only")
out <- do.call(rbind, res); print(out)
write.csv(out, "results/tcga/n_sig_summary.csv")
