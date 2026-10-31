# Volcano plot and 3-tool overlap plot for the positive control: Basal vs LumA (basal = case).
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
d <- read.csv("results/tcga/de_posctrl_basal_vs_lumA.csv")
tools <- c(DESeq2 = "q_deseq2", edgeR = "q_edger", `limma-voom` = "q_limmavoom")

# ---- Volcano (DESeq2 q-values) ----
d$q <- pmax(d$q_deseq2, 1e-300)
d$nlq <- -log10(d$q)
d$status <- ifelse(!is.na(d$q_deseq2) & d$q_deseq2 < 0.05 & abs(d$lfc) >= 1,
                   ifelse(d$lfc > 0, "Higher in basal", "Higher in LumA"), "Not significant")
markers <- c("KRT5", "KRT14", "KRT17", "FOXC1", "SOX10", "ESR1", "GATA3", "FOXA1", "PGR", "XBP1")
lab <- d[d$gene_name %in% markers & !is.na(d$q_deseq2), ]
cat("Known marker genes (positive lfc = higher in basal):\n")
print(data.frame(gene = lab$gene_name, lfc = round(lab$lfc, 2), q_deseq2 = signif(lab$q_deseq2, 2)), row.names = FALSE)
p <- ggplot(d[!is.na(d$nlq), ], aes(lfc, nlq, colour = status)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = c("Higher in basal" = "red", "Higher in LumA" = "blue", "Not significant" = "grey70")) +
  geom_text_repel(data = lab, aes(label = gene_name), colour = "black", size = 3, max.overlaps = Inf) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(x = "log2 fold change (as saved in the file), basal vs LumA", y = "-log10 q-value (DESeq2)",
       title = "Basal vs LumA: does the pipeline recover known biology?", colour = NULL) +
  theme_bw()
ggsave("results/figures/volcano_basal_vs_lumA.png", p, width = 8, height = 6, dpi = 150)

# ---- Overlap of genes called at q < 0.05 (UpSet-style, base R) ----
sig <- sapply(tools, function(n) !is.na(d[[n]]) & d[[n]] < 0.05)
cat("\nGenes called at q < 0.05 per tool:\n"); print(colSums(sig))
key <- apply(sig, 1, function(r) paste(as.integer(r), collapse = ""))
tab <- sort(table(key[key != "000"]), decreasing = TRUE)
cat("\nOverlap table (bits = DESeq2, edgeR, limma-voom):\n"); print(tab)
write.csv(data.frame(combo = names(tab), genes = as.vector(tab)), "results/tcga/posctrl_overlap_counts.csv", row.names = FALSE)
n <- length(tab); xl <- c(0.2, 1.2 * n)
png("results/figures/upset_basal_vs_lumA.png", width = 1400, height = 900, res = 150)
layout(matrix(1:2, nrow = 2), heights = c(3, 1.4))
par(mar = c(0.5, 7, 3, 1), xaxs = "i")
bp <- barplot(as.vector(tab), names.arg = NA, xlim = xl, ylab = "Genes", main = "Genes called by each combination of tools (q < 0.05)")
text(bp, as.vector(tab), labels = as.vector(tab), pos = 3, cex = 0.7, xpd = TRUE)
par(mar = c(3, 7, 0.5, 1), xaxs = "i")
plot(NA, xlim = xl, ylim = c(0.5, 3.5), axes = FALSE, xlab = "", ylab = "")
axis(2, at = 3:1, labels = names(tools), las = 1, tick = FALSE)
for (j in seq_len(n)) {
  bits <- as.integer(strsplit(names(tab)[j], "")[[1]])
  points(rep(bp[j], 3), 3:1, pch = 19, cex = 1.6, col = ifelse(bits == 1, "black", "grey88"))
  on <- (3:1)[bits == 1]; if (length(on) > 1) lines(rep(bp[j], 2), range(on), lwd = 2)
}
dev.off()
cat("\nsaved volcano_basal_vs_lumA.png and upset_basal_vs_lumA.png in results/figures\n")
