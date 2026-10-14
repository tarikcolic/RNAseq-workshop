###############################################################################
# Parathyroid adenoma RNA-seq - differential expression with DESeq2
# (translated from the cichlid pipeline)
#
# Design: patient x treatment x time (23 libraries, 4 patients)
#
# Analyses:
#   1a DPN vs Control, design ~ patient + treatment              (model WITHOUT time)
#   1b DPN vs Control, design ~ patient + time + treatment       (model WITH time)
#       -> compare the number of significant genes, decide which to publish
#   2  OHT vs Control, design ~ patient + time + treatment       (expected: ~nothing)
#
# Set the working directory to the folder containing the three data files.
###############################################################################

# =============================================================================
# SETUP
# =============================================================================
BiocManager::install(c("DESeq2", "apeglm", "org.Hs.eg.db", "AnnotationDbi"))
library(DESeq2)
library(apeglm)
library(org.Hs.eg.db)
library(AnnotationDbi)

out_dir <- "analyses/results_DESeq2_parathyroid"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# READ THE DATA
# =============================================================================
# Both files are comma separated -> read.csv(). First column = row names.
counts  <- read.csv("Data/parathyroid_counts.csv",  row.names = 1, check.names = FALSE)
samples <- read.csv("Data/parathyroid_samples.csv", row.names = 1)

dim(counts)            # expect 63193 x 23
head(samples)

#counts <- as.matrix(counts) # did not need

# patient is stored as an INTEGER: it must be a factor, otherwise DESeq2
# would treat it as a continuous co-variate (one slope instead of 4 levels).
# First level = reference. Control and 24h are the references.
samples$patient   <- factor(samples$patient)
samples$treatment <- factor(samples$treatment, levels = c("Control", "DPN", "OHT"))
samples$time      <- factor(samples$time,      levels = c("24h", "48h"))

# Group sizes. Note patient 4 has no Control_24h -> design is unbalanced.
table(samples$treatment, samples$time)
table(samples$patient, samples$treatment)

# Counts columns and sample rows MUST match. Reorder defensively, then check.
counts <- counts[, rownames(samples)]
stopifnot(all(colnames(counts) == rownames(samples)))

# Gene annotation: keep ENSEMBL IDs as row names (symbols are not unique and
# some are NA), and store the symbol as an extra column to add to the output.
symbols <- mapIds(org.Hs.eg.db, keys = rownames(counts),
                  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

# Helper: attach symbol and save
save_res <- function(res, file) {
  df <- as.data.frame(res)
  df$symbol <- symbols[rownames(df)]
  write.csv(df, file.path(out_dir, file))
  invisible(df)
}

# Filter rule (same logic as the cichlid script): >= 10 counts in at least as
# many samples as the smallest group. Smallest treatment group = Control, n = 7.
min_n <- min(table(samples$treatment))   # 7

# =============================================================================
# ANALYSIS 1a - DPN vs Control, WITHOUT time:  ~ patient + treatment
# =============================================================================
dds1a <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = samples,
                                design = ~ patient + treatment)

dds1a <- dds1a[rowSums(counts(dds1a) >= 10) >= min_n, ]
dds1a <- DESeq(dds1a)
resultsNames(dds1a)

res1a <- results(dds1a, contrast = c("treatment", "DPN", "Control"), alpha = 0.05)
summary(res1a)
n1a <- sum(res1a$padj < 0.05, na.rm = TRUE)
n1a

# =============================================================================
# ANALYSIS 1b - DPN vs Control, WITH time:  ~ patient + time + treatment
# The variable of interest goes LAST.
# =============================================================================
dds1b <- DESeqDataSetFromMatrix(countData = counts, colData = samples,
                                design = ~ patient + time + treatment)
dds1b <- dds1b[rowSums(counts(dds1b) >= 10) >= min_n, ]
dds1b <- DESeq(dds1b)
resultsNames(dds1b)

# =============================================================================
# QUALITY CONTROL (on the full model, dds1b) - do this before results()
# =============================================================================

# ---- QC 1. Size factors ------------------------------------------------------
sizeFactors(dds1b)
summary(sizeFactors(dds1b))
max(sizeFactors(dds1b)) / min(sizeFactors(dds1b))
barplot(sizeFactors(dds1b), las = 2, cex.names = 0.5,
        ylab = "size factor", main = "Size factors")
abline(h = 1, col = "red")

plot(colSums(counts(dds1b)) / 1e6, sizeFactors(dds1b),
     xlab = "library size (million counts)", ylab = "size factor",
     pch = 16, main = "Size factor vs library size")

boxplot(log2(counts(dds1b, normalized = TRUE) + 1), las = 2, cex.axis = 0.5,
        ylab = "log2(normalised count + 1)", main = "Normalised counts")

# ---- QC 2. Dispersion --------------------------------------------------------
plotDispEsts(dds1b)
summary(dispersions(dds1b))

# ---- QC 3. Mean-variance before / after VST ---------------------------------
norm_counts <- counts(dds1b, normalized = TRUE)
plot(rowMeans(norm_counts), apply(norm_counts, 1, var), log = "xy",
     pch = 16, cex = 0.3, xlab = "mean normalised count", ylab = "variance",
     main = "Raw counts: variance depends on the mean")
abline(0, 1, col = "red")

vsd <- vst(dds1b, blind = FALSE)
vst_mat <- assay(vsd)
plot(rowMeans(vst_mat), apply(vst_mat, 1, sd), pch = 16, cex = 0.3,
     xlab = "mean of VST values", ylab = "standard deviation",
     main = "After VST: spread is roughly flat")

# ---- QC 4. Cook's distance ---------------------------------------------------
# Only 2-4 replicates per cell, so DESeq2 will not auto-replace outliers
# (needs >= 7 per group): just note any sample with a much higher box.
boxplot(log10(assays(dds1b)[["cooks"]]), range = 0, las = 2, cex.axis = 0.5,
        ylab = "log10 Cook's distance", main = "Cook's distance per sample")

# ---- QC 5. Sparsity ----------------------------------------------------------
plotSparsity(dds1b)

# ---- QC 6. PCA: does the biggest axis match the design? ---------------------
# Expectation: PC1 = time (the large time-in-culture effect), PC2 = patient,
# treatment barely visible. That is exactly why time belongs in the model.
plotPCA(vsd, intgroup = c("time", "treatment"))
plotPCA(vsd, intgroup = c("patient", "time"))

sample_dist <- dist(t(vst_mat))
plot(hclust(sample_dist),
     labels = paste(samples$patient, samples$treatment, samples$time),
     cex = 0.6, main = "Sample clustering (VST, Euclidean)")

# =============================================================================
# ANALYSIS 1b - results
# =============================================================================
res1b <- results(dds1b, contrast = c("treatment", "DPN", "Control"), alpha = 0.05)
summary(res1b)
n1b <- sum(res1b$padj < 0.05, na.rm = TRUE)
n1b

# Compare the two models: same reads, different number of hits
c(without_time = n1a, with_time = n1b, fold_change = n1b / n1a)

# Shrunken fold changes (apeglm needs a coefficient; Control is the reference
# so the coefficient is treatment_DPN_vs_Control)
res1b_shrunk <- lfcShrink(dds1b, coef = "treatment_DPN_vs_Control", type = "apeglm")
res1b_shrunk <- res1b_shrunk[order(res1b_shrunk$padj), ]
head(res1b_shrunk)

# ---- QC 7. p-value histograms (raw p-values) --------------------------------
par(mfrow = c(1, 2))
hist(res1a$pvalue, breaks = 50, col = "grey", xlab = "raw p-value",
     main = "DPN vs Control: ~ patient + treatment")
hist(res1b$pvalue, breaks = 50, col = "grey", xlab = "raw p-value",
     main = "DPN vs Control: ~ patient + time + treatment")
par(mfrow = c(1, 1))

save_res(res1a,        "01a_DPN_vs_Control_no_time.csv")
save_res(res1b_shrunk, "01b_DPN_vs_Control_with_time.csv")

plotMA(res1b_shrunk, ylim = c(-5, 5))
plotCounts(dds1b, gene = rownames(res1b_shrunk)[1], intgroup = c("treatment", "time"))

# =============================================================================
# ANALYSIS 2 - OHT vs Control (same full model as 1b; an honest negative result)
# =============================================================================
resultsNames(dds1b)

res2 <- results(dds1b, contrast = c("treatment", "OHT", "Control"), alpha = 0.05)
summary(res2)
n2 <- sum(res2$padj < 0.05, na.rm = TRUE)
n2

hist(res2$pvalue, breaks = 50, col = "grey", xlab = "raw p-value",
     main = "OHT vs Control: p-value histogram")
# Expect an approximately flat histogram with little/no spike at 0.

res2_shrunk <- lfcShrink(dds1b, coef = "treatment_OHT_vs_Control", type = "apeglm")
res2_shrunk <- res2_shrunk[order(res2_shrunk$padj), ]
head(res2_shrunk)

save_res(res2_shrunk, "02_OHT_vs_Control.csv")

plotMA(res2_shrunk, ylim = c(-5, 5))

# =============================================================================
# OVERVIEW PCA
# =============================================================================
vsd_blind <- vst(dds1b, blind = TRUE)
png(file.path(out_dir, "00_PCA.png"), width = 1800, height = 1300, res = 240)
print(plotPCA(vsd_blind, intgroup = c("treatment", "time")))
dev.off()

# =============================================================================
# SUMMARY
# =============================================================================
c(DPN_vs_Control_no_time   = n1a,   # 1a
  DPN_vs_Control_with_time = n1b,   # 1b
  OHT_vs_Control           = n2)    # 2

# sessionInfo()
