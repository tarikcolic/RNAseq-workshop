###############################################################################
# Midas cichlid RNA-seq - differential expression with DESeq2
#
# Four analyses:
#   1  species:  cit vs ast
#   2  age:      d118 vs d7, separately within each species
#   3  interaction: does the age effect differ between species?
#   4  interaction: does the light effect differ between species? (d118 only)
#
# Run it line by line with Ctrl+Enter (Windows) or Cmd+Enter (Mac).
# Open 24_Sarajevo_Summer_School.Rproj first, so the file paths below work.
###############################################################################


# =============================================================================
# SETUP
# =============================================================================
BiocManager::install("apeglm", type = "binary")
library(DESeq2)
library(apeglm)

# Make a folder for the output
dir.create(
  "/Users/taric/Documents/GitHub/RNAseq-workshop/analyses/results_DESeq2", 
  recursive = TRUE, 
  showWarnings = FALSE
  )


# =============================================================================
# READ THE DATA
# =============================================================================

# The count file is TAB separated, so we use read.delim().
# row.names = 1 makes the first column (gene_id) into the row names.
counts <- read.delim(
  "/Users/taric/Documents/GitHub/RNAseq-workshop/Data/cichlid_counts.tsv", 
  row.names = 1)

# The sample sheet is COMMA separated, so we use read.csv().
samples <- read.csv(
  "/Users/taric/Documents/GitHub/RNAseq-workshop/Data/cichlid_samples.csv", 
  row.names = 1)
genes <- read.delim(
  "/Users/taric/Documents/GitHub/RNAseq-workshop/Day2/cichlid_geneNames.txt", 
  row.names = 1)

View(counts)
View(samples)
View(genes)

# Replace the obscure gene names
row.names(counts) <- genes$GeneName
# Look at them
dim(counts)
counts[1:5, 1:5]
head(samples)

# DESeq2 needs a matrix of whole numbers
counts <- as.matrix(counts)

# Turn the experimental variables into factors.
# The FIRST level is the reference - everything is compared against it.
samples$species <- factor(samples$species, levels = c("ast", "cit"))
samples$age     <- factor(samples$age,     levels = c("d7", "d118"))
samples$light   <- factor(samples$light,   levels = c("white", "red"))

# How many samples in each group?
table(samples$species, samples$age, samples$light)

# The columns of the counts and the rows of the sample sheet MUST be in the
# same order. If this prints FALSE, stop and fix it before going on.
all(colnames(counts) == rownames(samples))



# =============================================================================
# ANALYSIS 1 - SPECIES:  cit vs ast
#
# design = ~ light + age + species
# The variable of interest goes LAST. The others are there to account for
# the fact that samples also differ in light and age.
# =============================================================================

dds1 <- DESeqDataSetFromMatrix(countData = counts,
            colData   = samples,
            design    = ~ light + age + species)

str(dds1)


# Remove genes with almost no reads. The rule: keep a gene if it has at
# least 10 counts in at least as many samples as the smallest group.
# The smallest species group has 22 samples.
keep1 <- rowSums(counts(dds1) >= 10) >= 22
table(keep1)
dds1 <- dds1[keep1, ]

# Fit the model. This does the normalization, dispersion estimation and
# testing in one step.
dds1 <- DESeq(dds1)

# =============================================================================
# QUALITY CONTROL OF THE FITTED OBJECT
#
# Do this EVERY TIME, after DESeq() and BEFORE you look at results().
# If the model does not fit, the p-values are meaningless no matter how
# small they are.
#
# DESeq2 has exactly five plotting functions:
#     plotDispEsts   dispersion estimates      <- the main fit check
#     plotMA         fold change vs expression
#     plotPCA        samples in 2D
#     plotCounts     one gene at a time
#     plotSparsity   are counts concentrated in a few samples?
# Everything else below uses base R on numbers DESeq2 stores in the object.
# =============================================================================

# ---- QC 1. Size factors: did normalization work? ----------------------------
# One number per sample. DESeq2 divides each column by its size factor.
sizeFactors(dds1)

# They should sit near 1. A big spread means the libraries differ a lot in
# depth - not fatal, but worth knowing, because the shallowest samples carry
# the least information.
# In THIS dataset they run 0.43 to 2.14, a 4.9-fold spread, which is wide.
# That follows from the shallow sequencing (median ~1.3 M counts per sample).
summary(sizeFactors(dds1))

max(sizeFactors(dds1)) / min(sizeFactors(dds1))

barplot(sizeFactors(dds1), las = 2, cex.names = 0.5,
     ylab = "size factor", main = "Size factors")

abline(h = 1, col = "red")

# Size factors should track the raw library sizes. If a point is far off the
# line, that sample has an unusual count distribution, not just more reads.
plot(colSums(counts(dds1)) / 1e6, sizeFactors(dds1),
   xlab = "library size (million counts)", ylab = "size factor",
   pch = 16, main = "Size factor vs library size")

# Box plot of normalized counts on a log scale. After normalization the
# boxes should line up. normalized = TRUE is what applies the size factors.
boxplot(log2(counts(dds1, normalized = TRUE) + 1), las = 2, cex.axis = 0.5,
    ylab = "log2(normalised count + 1)", main = "Normalised counts")


# ---- QC 2. plotDispEsts: The dispersion plot --------------------------------
plotDispEsts(dds1)

# How to read it. Each dot is one gene, x is mean normalised count,
# y is dispersion (how much more variable the gene is than Poisson noise).
#
#   BLACK dots  the raw per-gene estimate, from that gene's data alone
#   RED line    the fitted trend: typical dispersion at that expression level
#   BLUE dots   the final shrunken estimate, pulled toward the red line
#   BLUE circles genes NOT shrunk, flagged as dispersion outliers
#
# What a GOOD plot looks like:
#   - the red line falls as expression rises, then flattens
#   - black dots scatter around the red line
#   - blue dots sit close to the line but keep some spread
#
# What a PROBLEM looks like:
#   - black dots mostly ABOVE the line -> unmodelled variation. You are
#     missing a covariate (a batch, a sex, a plate) from the design.
#   - the red line is flat and high -> very noisy data
#   - almost no spread in the blue dots -> over-shrunk, too few replicates
#
# The same numbers, if you want them rather than the picture:
summary(dispersions(dds1))
dispersionFunction(dds1)


# ---- QC 3. Variance estimation: mean-variance before and after -------------
# This is the problem the dispersion model exists to solve: in raw counts
# the variance grows with the mean, so highly expressed genes dominate any
# distance or clustering calculation.

norm_counts <- counts(dds1, normalized = TRUE)
gene_mean <- rowMeans(norm_counts)
gene_var  <- apply(norm_counts, 1, var)

# On raw normalised counts: variance tracks the mean (slope about 1 on a
# log-log plot). This is the mean-variance dependence.
plot(gene_mean, gene_var, log = "xy", pch = 16, cex = 0.3,
    xlab = "mean normalised count", ylab = "variance",
    main = "Raw counts: variance depends on the mean")

abline(0, 1, col = "red")       # the Poisson expectation, variance = mean

# vst() removes that dependence. After transforming, the spread should be
# roughly flat across expression levels.
vsd1 <- vst(dds1, blind = FALSE)
vst_mat <- assay(vsd1)
plot(rowMeans(vst_mat), apply(vst_mat, 1, sd), pch = 16, cex = 0.3,
   xlab = "mean of VST values", ylab = "standard deviation",
   main = "After VST: spread is roughly flat")

# Use vsd for PCA, clustering and heatmaps. Use the raw counts in dds for
# results(). Never swap them around.


# ---- QC 4. Cook's distance: is one sample driving a gene? -------------------
# DESeq2 stores a Cook's distance per gene per sample. A high value means
# that one sample alone is pulling that gene's fold change.
boxplot(log10(assays(dds1)[["cooks"]]), range = 0, las = 2, cex.axis = 0.5,
     ylab = "log10 Cook's distance", main = "Cook's distance per sample")

# One box much higher than the rest = a sample with many influential genes,
# i.e. a candidate outlier. With 7 or more replicates per group DESeq2
# automatically replaces such counts, so you usually just note it here.


# ---- QC 5. plotSparsity: are counts concentrated in few samples? -----------
plotSparsity(dds1)

# x is the total count for a gene, y is the fraction of that total coming
# from its single largest sample. Points near y = 1 on the right are genes
# whose entire signal is one sample - often technical. A cloud pushed
# towards the top means the data may not suit a count-based GLM.


# ---- QC 6. plotPCA: do the samples group as the design says? ---------------
plotPCA(vsd1, intgroup = c("species", "age"))

# Then ask: does the biggest axis of variation match the biggest factor in
# my design? Here PC1 should be age. If instead PC1 separated samples by
# sequencing run, that would need to be in the design formula.

# Sample-to-sample distances, using base R only
sample_dist <- dist(t(assay(vsd1)))

plot(hclust(sample_dist), labels = paste(samples$species, samples$age),
   cex = 0.6, main = "Sample clustering (VST, Euclidean)")


## Save the QC figures to one PDF, so you have a record
#plotDispEsts(dds1)
#plotSparsity(dds1)
#plotPCA(vsd1, intgroup = c("species", "age"))
#barplot(sizeFactors(dds1), las = 2, cex.names = 0.5, ylab = "size factor")
#boxplot(log2(counts(dds1, normalized = TRUE) + 1), las = 2, cex.axis = 0.5,
#     main = "Normalised counts")
#boxplot(log10(assays(dds1)[["cooks"]]), range = 0, las = 2, cex.axis = 0.5,
#     main = "Cook's distance")




# What comparisons are available?
resultsNames(dds1)

# Get the results for species
res1 <- results(dds1, name = "species_cit_vs_ast", alpha = 0.05)
summary(res1)

# How many genes are significant?
sum(res1$padj < 0.05, na.rm = TRUE)

# Shrink the log2 fold changes. This makes the fold changes of lowly
# expressed genes less extreme, which is better for ranking and plotting.
res1_shrunk <- lfcShrink(dds1, coef = "species_cit_vs_ast", type = "apeglm")

# Sort by adjusted p-value and look at the top genes
res1_shrunk <- res1_shrunk[order(res1_shrunk$padj), ]
head(res1_shrunk)


# ---- QC 7. The p-value histogram -------------------------------------------
# The most informative single plot in the whole analysis. Use the RAW
# p-values, not padj.
hist(res1$pvalue, breaks = 50, col = "grey",
   xlab = "raw p-value", main = "p-value histogram")

# What you want: a flat carpet across the middle, with a spike at zero.
# The carpet is the genes with no effect; the spike is your signal.
#
# Shapes that mean trouble:
#   hump in the MIDDLE (around 0.5)  -> the model is wrong, usually a
#                                       missing covariate
#   spike at 1 as well as 0          -> over-dispersed, or filtering issue
#   completely flat, no spike at 0   -> no detectable effect (a real result)
#   spiky / comb-like                -> very low counts, discreteness


# Save
write.csv(as.data.frame(res1_shrunk), "analyses/results_DESeq2/01_species_cit_vs_ast.csv")

# Two quick diagnostic plots
plotMA(res1_shrunk, ylim = c(-5, 5))
plotCounts(dds1, gene = rownames(res1_shrunk)[1], intgroup = "species")



# =============================================================================
# ANALYSIS 2a - AGE WITHIN ast:  d118 vs d7
#
# We keep only the ast samples, then compare ages inside that group.
# =============================================================================

# Which samples are ast?
samples_ast <- samples[samples$species == "ast", ]
counts_ast  <- counts[, rownames(samples_ast)]

# droplevels() removes the now-unused "cit" level from the species factor
samples_ast <- droplevels(samples_ast)
table(samples_ast$age, samples_ast$light)

dds2a <- DESeqDataSetFromMatrix(countData = counts_ast,
  colData   = samples_ast,
  design    = ~ light + age)

# Smallest age group within ast has 12 samples
keep2a <- rowSums(counts(dds2a) >= 10) >= 12
dds2a  <- dds2a[keep2a, ]

dds2a <- DESeq(dds2a)
resultsNames(dds2a)

res2a <- results(dds2a, name = "age_d118_vs_d7", alpha = 0.05)
summary(res2a)
sum(res2a$padj < 0.05, na.rm = TRUE)

res2a_shrunk <- lfcShrink(dds2a, coef = "age_d118_vs_d7", type = "apeglm")
res2a_shrunk <- res2a_shrunk[order(res2a_shrunk$padj), ]
head(res2a_shrunk)

write.csv(as.data.frame(res2a_shrunk),
                       "analyses/results_DESeq2/02a_age_d118_vs_d7_in_ast.csv")

plotMA(res2a_shrunk, ylim = c(-5, 5))
plotCounts(dds2a, gene = rownames(res2a_shrunk)[1], intgroup = "age")



# =============================================================================
# ANALYSIS 2b - AGE WITHIN cit:  d118 vs d7
#
# Exactly the same steps as 2a, with cit instead of ast.
# =============================================================================

samples_cit <- samples[samples$species == "cit", ]
counts_cit  <- counts[, rownames(samples_cit)]
samples_cit <- droplevels(samples_cit)
table(samples_cit$age, samples_cit$light)

dds2b <- DESeqDataSetFromMatrix(countData = counts_cit,
            colData   = samples_cit,
            design    = ~ light + age)

# Smallest age group within cit has 10 samples
keep2b <- rowSums(counts(dds2b) >= 10) >= 10
dds2b  <- dds2b[keep2b, ]

dds2b <- DESeq(dds2b)
resultsNames(dds2b)

res2b <- results(dds2b, name = "age_d118_vs_d7", alpha = 0.05)
summary(res2b)
sum(res2b$padj < 0.05, na.rm = TRUE)

res2b_shrunk <- lfcShrink(dds2b, coef = "age_d118_vs_d7", type = "apeglm")
res2b_shrunk <- res2b_shrunk[order(res2b_shrunk$padj), ]
head(res2b_shrunk)

write.csv(as.data.frame(res2b_shrunk),
     "analyses/results_DESeq2/02b_age_d118_vs_d7_in_cit.csv")

plotMA(res2b_shrunk, ylim = c(-5, 5))
plotCounts(dds2b, gene = rownames(res2b_shrunk)[1], intgroup = "age")


# =============================================================================
# ANALYSIS 3 - INTERACTION:  species x age
#
# The question is NOT "does age matter" and NOT "do species differ".
# It is: IS THE AGE EFFECT THE SAME IN BOTH SPECIES?
#
# We add the term species:age to the design. Its coefficient is the
# DIFFERENCE between the age effect in cit and the age effect in ast.
# =============================================================================

dds3 <- DESeqDataSetFromMatrix(countData = counts,
                        colData   = samples,
                        design    = ~ light + species + age + species:age)

keep3 <- rowSums(counts(dds3) >= 10) >= 22
dds3  <- dds3[keep3, ]

dds3 <- DESeq(dds3)

# The interaction is the last name in this list
resultsNames(dds3)

res3 <- results(dds3, name = "speciescit.aged118", alpha = 0.05)
summary(res3)
sum(res3$padj < 0.05, na.rm = TRUE)

res3 <- res3[order(res3$padj), ]
head(res3)

write.csv(as.data.frame(res3),
        "analyses/results_DESeq2/03_interaction_species_x_age.csv")

# Look at the top interaction gene. The two lines should NOT be parallel -
# that is what an interaction means.
plotCounts(dds3, gene = rownames(res3)[1], intgroup = c("age", "species"))



# =============================================================================
# ANALYSIS 4 - INTERACTION:  species x light, using ONLY the d118 samples
#
# Question: do the two species respond differently to red vs white light?
# =============================================================================

samples_d118 <- samples[samples$age == "d118", ]
counts_d118  <- counts[, rownames(samples_d118)]
samples_d118 <- droplevels(samples_d118)

# Note the group sizes: cit/white has only 4 samples, the rest have 6.
table(samples_d118$species, samples_d118$light)

dds4 <- DESeqDataSetFromMatrix(countData = counts_d118,
                  colData   = samples_d118,
                  design    = ~ species + light + species:light)

# Smallest light group in this subset has 10 samples
keep4 <- rowSums(counts(dds4) >= 10) >= 10
dds4  <- dds4[keep4, ]

dds4 <- DESeq(dds4)
resultsNames(dds4)

res4 <- results(dds4, name = "speciescit.lightred", alpha = 0.05)
summary(res4)
sum(res4$padj < 0.05, na.rm = TRUE)

res4 <- res4[order(res4$padj), ]
head(res4)

write.csv(as.data.frame(res4),
      "analyses/results_DESeq2/04_interaction_species_x_light_d118.csv")

plotCounts(dds4, gene = rownames(res4)[1], intgroup = c("light", "species"))



# =============================================================================
# A PCA PLOT, TO SEE THE WHOLE EXPERIMENT AT ONCE
# =============================================================================

# vst() puts the counts on a log-like scale that is suitable for plotting.
# Never use it as input to results() - DESeq2 needs raw counts for that.
vsd <- vst(dds1, blind = TRUE)
plotPCA(vsd, intgroup = c("species", "age"))

png("analyses/results_DESeq2/00_PCA.png", width = 1800, height = 1300, res = 240)
print(plotPCA(vsd, intgroup = c("species", "age")))
dev.off()



# =============================================================================
# SUMMARY OF ALL FOUR ANALYSES
# =============================================================================

sum(res1$padj  < 0.05, na.rm = TRUE)   # 1  species
sum(res2a$padj < 0.05, na.rm = TRUE)   # 2a age within ast
sum(res2b$padj < 0.05, na.rm = TRUE)   # 2b age within cit
sum(res3$padj  < 0.05, na.rm = TRUE)   # 3  species x age
sum(res4$padj  < 0.05, na.rm = TRUE)   # 4  species x light, d118

#sessionInfo()