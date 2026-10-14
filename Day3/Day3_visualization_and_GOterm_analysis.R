###############################################################################
# Day 3 - PLOTS AND GO-TERM ANALYSIS, simple version
#
# Dataset: airway (human airway smooth muscle cells, dexamethasone vs control,
#          4 donors, 8 samples). ENSEMBL gene IDs, so org.Hs.eg.db works.
#
# What this script makes:
#   1  VST                        10 genes before and after transformation
#   2  PCA                        01_PCA.png
#   3  Sample distance heatmap    02_sample_distances.png
#   4  Volcano plot               03_volcano.png
#   5  Two heatmaps               04_heatmap_variable.png, 05_heatmap_significant.png
#   6  Gene ID translation        ENSEMBL -> SYMBOL -> ENTREZ
#   7  GO over-representation     06_GO_dotplot_up.png, 07_GO_barplot_down.png,
#                                 08_GO_emapplot.png, 09_GO_cnetplot.png
#   8  GSEA                       10_GSEA_dotplot.png, 11_GSEA_ridgeplot.png,
#                                 12_GSEA_runningscore.png
#   9  Composite figure           13_composite.png
#
# Run it line by line with Ctrl+Enter (Windows) or Cmd+Enter (Mac).
###############################################################################


# =============================================================================
# 0. SETUP
# =============================================================================
install.packages("ape")
install.packages("ggiraph")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db", type = "source")

options(timeout = 1200)

library(ape)
library(ggiraph)
library(clusterProfiler)

library(DESeq2)
library(ggplot2)
library(ggrepel)         # for non-overlapping gene labels on the volcano
library(pheatmap)        # for heatmaps
library(RColorBrewer)    # for colour palettes
library(clusterProfiler) # for the GO analysis
library(org.Hs.eg.db)    # the human annotation database
library(enrichplot)      # for the GO plots
library(patchwork)       # for sticking plots together

# ---- Read the data ----------------------------------------------------------
counts  <- read.csv("/Users/taric/Documents/GitHub/RNAseq-workshop/Data/airway_counts.csv", row.names = 1)
samples <- read.csv("/Users/taric/Documents/GitHub/RNAseq-workshop/Data/airway_samples.csv", row.names = 1)

# The counts file has a gene_name column as well as the samples - drop it,
# DESeq2 needs numbers only.
counts <- counts[, rownames(samples)]
counts <- as.matrix(counts)

# dex is the treatment. untrt is the reference level.
samples$dex  <- factor(samples$dex,  levels = c("untrt", "trt"))
samples$cell <- factor(samples$cell)

table(samples$cell, samples$dex)
all(colnames(counts) == rownames(samples))


# ---- Fit the model ----------------------------------------------------------
# cell is in the design because the same 4 donors appear in both groups.
# This is a paired design.
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = samples,
                              design    = ~ cell + dex)

keep <- rowSums(counts(dds) >= 10) >= 1
table(keep)
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
summary(res)

# Shrunken fold changes, for plotting and ranking
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")



# =============================================================================
# 1. VARIANCE-STABILISING TRANSFORMATION
#
# Counts have variance that grows with the mean. Every plot in this script
# that measures DISTANCE between samples (PCA, heatmaps, clustering) needs
# that removed first, or the loudest genes decide the answer.
# =============================================================================

vsd <- vst(dds, blind = FALSE)

# Compare the two scales on the same two samples
par(mfrow = c(1, 2))
plot(log2(counts(dds, normalized = TRUE)[, 1:2] + 1), pch = 16, cex = 0.3,
     main = "log2(count + 1)")
plot(assay(vsd)[, 1:2], pch = 16, cex = 0.3, main = "VST")
par(mfrow = c(1, 1))

# Note what changed: the fan shape at low counts is gone.



# =============================================================================
# 2. PCA
# =============================================================================

# The quick version, straight from DESeq2
plotPCA(vsd, intgroup = "dex")

# The version you can customise: ask for the data instead of the plot
pca_data <- plotPCA(vsd, intgroup = "dex", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
head(pca_data)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = dex, shape = cell)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("airway: PCA on VST data") +
  theme_bw()

p_pca
ggsave("analyses/results_airway/01_PCA.png", plot = p_pca,
       width = 6, height = 4.5, dpi = 300)

# Read it: PC1 should separate the donors (cell), PC2 the treatment, or the
# other way round. If treatment is not on either of the first two PCs, the
# treatment effect is smaller than the donor-to-donor variation.



# =============================================================================
# 3. SAMPLE-DISTANCE HEATMAP
#
# Same information as the PCA, shown as all pairwise distances.
# =============================================================================

# dist() works on ROWS, and we want distances between SAMPLES, so t()
sample_dist <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dist)

# Label the rows with something readable
rownames(dist_matrix) <- paste(vsd$dex, vsd$cell, sep = " - ")
colnames(dist_matrix) <- NULL

colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colours,
         main = "Sample-to-sample distances",
         filename = "/Users/taric/Documents/GitHub/RNAseq-workshop/Day3/02_sample_distances.png",
         width = 6.5, height = 5)

# pheatmap writes straight to the file when you give it filename=.
# Leave filename out to see it on screen instead.
pheatmap(dist_matrix, clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist, col = colours,
         main = "Sample-to-sample distances")



# =============================================================================
# 4. VOLCANO PLOT
# =============================================================================

# Turn the results into a plain data frame
res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)

# Drop genes with no p-value (filtered out by DESeq2)
res_df <- res_df[!is.na(res_df$padj), ]

# A column saying what each gene is, for the colours
res_df$direction <- "not significant"
res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "up"
res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "down"
table(res_df$direction)

# ENSG00000165995 means nothing to a reader, so get gene symbols for labels.
# mapIds returns one symbol per gene ID, or NA.
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys      = res_df$gene,
                        keytype   = "ENSEMBL",
                        column    = "SYMBOL",
                        multiVals = "first")

# How many genes got a name?
mean(!is.na(res_df$symbol))

# Where the symbol is missing, fall back to the ENSEMBL ID
res_df$label <- ifelse(is.na(res_df$symbol), res_df$gene, res_df$symbol)

# Label only the 15 most significant genes - more than that is unreadable
res_df <- res_df[order(res_df$padj), ]
to_label <- res_df[1:15, ]

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = direction), size = 1, alpha = 0.6) +
  scale_colour_manual(values = c("up" = "firebrick",
                                 "down" = "steelblue",
                                 "not significant" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = to_label, aes(label = label),
                  size = 3, fontface = "italic", max.overlaps = 20) +
  labs(x = "log2 fold change (trt / untrt)", y = "-log10 adjusted p-value",
       title = "airway: dexamethasone vs control") +
  theme_bw()

p_volcano
ggsave("analyses/results_airway/03_volcano.png", plot = p_volcano,
       width = 7, height = 5.5, dpi = 300)



# =============================================================================
# 5. TWO HEATMAPS, AND WHY ONLY ONE OF THEM IS EVIDENCE
# =============================================================================

# ---- 5a. The top 30 most VARIABLE genes ------------------------------------
# Chosen without looking at the groups. If these genes still separate
# treated from control, that is real structure in the data.
gene_var <- rowVars(assay(vsd))
top_variable <- head(order(gene_var, decreasing = TRUE), 30)

mat_var <- assay(vsd)[top_variable, ]
mat_var <- mat_var - rowMeans(mat_var)      # centre each gene

# Use symbols as row names where we have them
rownames(mat_var) <- mapIds(org.Hs.eg.db, keys = rownames(mat_var),
                            keytype = "ENSEMBL", column = "SYMBOL",
                            multiVals = "first")
rownames(mat_var)[is.na(rownames(mat_var))] <- "(no symbol)"

annot <- as.data.frame(colData(vsd)[, c("dex", "cell")])

pheatmap(mat_var, annotation_col = annot,
         main = "Top 30 most variable genes (groups NOT used to choose them)",
         width = 7, height = 7)


# ---- 5b. The top 30 most SIGNIFICANT genes ---------------------------------
# Chosen BECAUSE they differ between the groups. This heatmap will always
# look clean - that is circular, not evidence. Show it to illustrate the
# effect, never as proof the effect exists.
top_sig <- head(rownames(res_df), 30)

mat_sig <- assay(vsd)[top_sig, ]
mat_sig <- mat_sig - rowMeans(mat_sig)
rownames(mat_sig) <- res_df$label[1:30]

pheatmap(mat_sig, annotation_col = annot,
         main = "Top 30 significant genes (circular - for illustration only)",
         width = 7, height = 7)



# =============================================================================
# 6. GENE IDENTIFIER TRANSLATION
#
# clusterProfiler's GO functions want ENTREZ IDs. Our data has ENSEMBL IDs.
# bitr() = "biological ID translator".
# =============================================================================

id_map <- bitr(res_df$gene,
               fromType = "ENSEMBL",
               toType   = "ENTREZID",
               OrgDb    = org.Hs.eg.db)

head(id_map)
nrow(id_map)

# Not every gene translates. Always check how many you lost - if it is more
# than about a third, something is wrong with your ID type.
nrow(id_map) / nrow(res_df)

# Add the ENTREZ IDs onto the results table.
#
# CAREFUL HERE. The obvious thing to write is
#     res_df <- merge(res_df, id_map, by.x = "gene", by.y = "ENSEMBL")
# but merge() keeps only the rows present in BOTH tables, so it silently
# DELETES every gene that has no ENTREZ ID - about 12% of them. Your
# results table then shrinks, and any count you report from it afterwards
# is too low. No warning is printed.
#
# Use match() instead: it adds a column and leaves the rows alone,
# putting NA where there is no match.
res_df$ENTREZID <- id_map$ENTREZID[match(res_df$gene, id_map$ENSEMBL)]

# Check the damage the merge would have done
nrow(res_df)
sum(!is.na(res_df$ENTREZID))

res_df <- res_df[order(res_df$padj), ]



# =============================================================================
# 7. GO OVER-REPRESENTATION ANALYSIS
#
# The question: among my significant genes, is any GO term more common than
# you would expect by chance?
#
# THE UNIVERSE MATTERS. The comparison set must be the genes that COULD
# have been detected in this experiment, not every gene in the genome.
# Getting this wrong is the most common mistake in pathway analysis.
# =============================================================================

# The three gene lists we need
sig_up   <- res_df$ENTREZID[res_df$padj < 0.05 & res_df$log2FoldChange >  1]
sig_down <- res_df$ENTREZID[res_df$padj < 0.05 & res_df$log2FoldChange < -1]
universe <- res_df$ENTREZID        # everything we tested

length(sig_up); length(sig_down); length(universe)


# ---- Up-regulated genes ------------------------------------------------------
ego_up <- enrichGO(gene          = sig_up,
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",          # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)          # show gene symbols

head(as.data.frame(ego_up)[, c("Description", "GeneRatio", "p.adjust")], 10)
nrow(as.data.frame(ego_up))

# GO terms overlap heavily - many are parent and child of each other.
# simplify() collapses near-duplicates so the plot is readable.
ego_up_simple <- simplify(ego_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
nrow(as.data.frame(ego_up_simple))

# GO term names are long. label_format sets how many characters before the
# label wraps to a new line - the default of 30 makes 3-line labels that
# collide with each other. 45 keeps them to two lines.
p_dot_up <- dotplot(ego_up_simple, showCategory = 12, label_format = 45) +
  ggtitle("GO Biological Process: UP in dexamethasone")

p_dot_up
ggsave("analyses/results_airway/06_GO_dotplot_up.png", plot = p_dot_up,
       width = 8.5, height = 7, dpi = 300)


# ---- Downregulated genes ---------------------------------------------------
ego_down <- enrichGO(gene = sig_down, universe = universe,
                     OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2, readable = TRUE)

ego_down_simple <- simplify(ego_down, cutoff = 0.7, by = "p.adjust",
                            select_fun = min)

p_bar_down <- barplot(ego_down_simple, showCategory = 12, label_format = 45) +
  ggtitle("GO Biological Process: DOWN in dexamethasone")

p_bar_down
ggsave("analyses/results_airway/07_GO_barplot_down.png", plot = p_bar_down,
       width = 8, height = 6.5, dpi = 300)


# ---- Enrichment map: which terms share genes? ------------------------------
# pairwise_termsim() must be run first - emapplot needs the similarities.
ego_up_sim <- pairwise_termsim(ego_up_simple)

p_emap <- emapplot(ego_up_sim, showCategory = 20) +
  ggtitle("GO terms that share genes (UP)")

p_emap
ggsave("analyses/results_airway/08_GO_emapplot.png", plot = p_emap,
       width = 8, height = 7, dpi = 300)


# ---- Gene-concept network: which genes drive which term? -------------------
# This is the plot that turns a term list back into biology.
#
# We colour the gene nodes by fold change. Because we used readable = TRUE,
# the object holds gene SYMBOLS, so the vector we pass must be named with
# symbols too.
fc_by_symbol <- res_df$log2FoldChange
names(fc_by_symbol) <- res_df$symbol
fc_by_symbol <- fc_by_symbol[!is.na(names(fc_by_symbol))]

p_cnet <- cnetplot(ego_up_simple, showCategory = 5,
                   foldChange = fc_by_symbol)

p_cnet
ggsave("analyses/results_airway/09_GO_cnetplot.png", plot = p_cnet,
       width = 9, height = 8, dpi = 300)


## ---- The universe mistake, demonstrated ------------------------------------
## Run the same thing with no universe. clusterProfiler then compares against
## all annotated genes in the genome, including ones not expressed in lung.
#ego_wrong <- enrichGO(gene = sig_up, universe = NULL,
#                      OrgDb = org.Hs.eg.db, ont = "BP",
#                      pvalueCutoff = 0.05, readable = TRUE)
#
## Compare the counts. The wrong version is usually LONGER and more exciting.
#nrow(as.data.frame(ego_up))
#nrow(as.data.frame(ego_wrong))



# =============================================================================
# 8. GENE SET ENRICHMENT ANALYSIS (GSEA)
#
# Over-representation needs a cutoff, which throws away everything just
# below it. GSEA uses the WHOLE ranked list instead, so a pathway where
# every gene shifts slightly can still be detected.
# =============================================================================

# Build the ranked vector: named numbers, sorted high to low.
# We rank by shrunken log2 fold change.
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$ENTREZID

# Remove duplicates and NAs, then sort - gseGO requires decreasing order
gene_list <- gene_list[!is.na(names(gene_list)) & !is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

head(gene_list)
length(gene_list)

# gseGO uses the OrgDb directly, so no download is needed.
set.seed(42)    # GSEA uses permutations, so fix the seed to get the same
                # answer twice
gsea <- gseGO(geneList     = gene_list,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              keyType      = "ENTREZID",
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

gsea_df <- as.data.frame(gsea)
nrow(gsea_df)
head(gsea_df[, c("Description", "setSize", "NES", "p.adjust")], 10)

# NES is the normalised enrichment score: positive = the set sits at the top
# of the ranking (up with treatment), negative = at the bottom.

# Dotplot split by direction
p_gsea_dot <- dotplot(gsea, showCategory = 10, split = ".sign", label_format = 40) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA: GO Biological Process")

p_gsea_dot
ggsave("analyses/results_airway/10_GSEA_dotplot.png", plot = p_gsea_dot,
       width = 9, height = 7, dpi = 300)

# Ridgeplot: the distribution of fold changes inside each set
p_ridge <- ridgeplot(gsea, showCategory = 15, label_format = 45) +
  ggtitle("Fold-change distribution within each gene set")

p_ridge
ggsave("analyses/results_airway/11_GSEA_ridgeplot.png", plot = p_ridge,
       width = 8, height = 7, dpi = 300)




# =============================================================================
# 9. SAVE THE TABLES
# =============================================================================

write.csv(res_df, "analyses/results_airway/DE_results.csv", row.names = FALSE)
write.csv(as.data.frame(ego_up_simple),
          "analyses/results_airway/GO_up.csv", row.names = FALSE)
write.csv(as.data.frame(ego_down_simple),
          "analyses/results_airway/GO_down.csv", row.names = FALSE)
write.csv(gsea_df, "analyses/results_airway/GSEA.csv", row.names = FALSE)


# =============================================================================
# 11. THE RESULT IN WORDS
# =============================================================================

sum(res_df$padj < 0.05)                                   # significant genes
sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1)  # and 2-fold or more
nrow(as.data.frame(ego_up_simple))                        # GO terms up
nrow(as.data.frame(ego_down_simple))                      # GO terms down
nrow(gsea_df)                                             # GSEA sets

sessionInfo()
