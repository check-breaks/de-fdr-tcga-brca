# PCA / QC on TCGA-BRCA counts, coloured by PAM50 subtype and by tissue source site.
suppressMessages({library(SummarizedExperiment); library(edgeR)})
se <- readRDS("data/tcga/brca_se.rds")
cd <- as.data.frame(colData(se))
counts <- assay(se, "unstranded")

# 1. keep genes with CPM > 1 in at least 20 tumors (drops mostly-zero genes)
keep <- rowSums(cpm(counts) > 1) >= 20
cat("genes kept:", sum(keep), "of", nrow(counts), "\n")

# 2. log2 counts-per-million: puts samples on a comparable, roughly symmetric scale
lc <- cpm(counts[keep, ], log = TRUE, prior.count = 1)

# 3. PCA on the 500 most variable genes
top <- order(apply(lc, 1, var), decreasing = TRUE)[1:500]
pc <- prcomp(t(lc[top, ]), scale. = TRUE)
pve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)

subtype <- as.character(cd$BRCA_Subtype_PAM50); subtype[is.na(subtype)] <- "NA"
site <- sub("^TCGA-([A-Z0-9]+)-.*$", "\\1", colnames(se))

# 4. How much of each PC does subtype / site explain? (R-squared from a simple linear model)
r2 <- function(f) sapply(1:5, function(i) round(summary(lm(pc$x[, i] ~ f))$r.squared, 2))
tab <- data.frame(PC = paste0("PC", 1:5), pct_variance = pve[1:5],
                  R2_subtype = r2(factor(subtype)), R2_site = r2(factor(site)))
print(tab)
write.csv(tab, "results/tcga/pca_r2_summary.csv", row.names = FALSE)

# 5. Figure: left = subtype (BRCA1 mutants ringed in black), right = top 8 sites
cols <- c(Basal = "red", Her2 = "purple", LumA = "blue", LumB = "skyblue", Normal = "darkgreen", "NA" = "grey70")
topsite <- names(sort(table(site), decreasing = TRUE))[1:8]
scol <- setNames(rainbow(8), topsite)
png("results/figures/pca_subtype_site.png", width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(pc$x[, 1], pc$x[, 2], pch = 19, col = cols[subtype], cex = 0.7,
     xlab = paste0("PC1 (", pve[1], "%)"), ylab = paste0("PC2 (", pve[2], "%)"),
     main = "Coloured by PAM50 subtype")
m <- which(cd$BRCA1_mut %in% TRUE)
points(pc$x[m, 1], pc$x[m, 2], pch = 1, cex = 1.4, lwd = 2)
legend("topright", legend = names(cols), col = cols, pch = 19, cex = 0.7)
plot(pc$x[, 1], pc$x[, 2], pch = 19, cex = 0.7,
     col = ifelse(site %in% topsite, scol[site], "grey80"),
     xlab = paste0("PC1 (", pve[1], "%)"), ylab = paste0("PC2 (", pve[2], "%)"),
     main = "Coloured by tissue source site (top 8)")
legend("topright", legend = topsite, col = scol, pch = 19, cex = 0.7)
dev.off()
cat("saved results/figures/pca_subtype_site.png\n")
