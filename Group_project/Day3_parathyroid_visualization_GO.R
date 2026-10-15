###############################################################################
# Day 3 for the parathyroid dataset - PLOTS AND GO-TERM ANALYSIS
# (adapted from Day3_visualization_and_GOterm_analysis.R, which used airway)
#
# Dataset: primary human parathyroid adenoma cells, 4 patients, DPN / OHT /
#          Control, 24 h and 48 h, 23 libraries. ENSEMBL IDs -> org.Hs.eg.db.
# Model:   ~ patient + time + treatment   (the model we decided to publish)
#
# What this script makes (all files go to out_dir):
#   1  Model fit + result tables   DPN vs Control and OHT vs Control
#   2  PCA                         01_PCA_raw_vs_cleaned.png
#   3  Sample distance heatmap     02_sample_distances.png
#   4  Volcano plots               03_volcano_DPN.png, 03_volcano_OHT.png
#   5  MA plots                    04_MA_DPN_OHT.png
#   6  Heatmaps                    05_heatmap_variable.png, 06_heatmap_DPN_significant.png
#   7  Single-gene plots           07_top_genes_DPN.png  (all 6 groups labelled)
#   8  Gene ID translation         ENSEMBL -> ENTREZ
#   9  GO over-representation      08_GO_dotplot_*.png, 09_GO_emapplot.png, 10_GO_cnetplot.png
#  10  GSEA (DPN and OHT)          11_GSEA_*.png
#  11  Composite figure            12_composite.png
#  12  Tables                      *.csv
#
# Run it line by line with Ctrl+Enter (Windows) or Cmd+Enter (Mac).
#
# THREE THINGS THAT DIFFER FROM THE AIRWAY TEMPLATE (and why)
#  a) No |log2FC| > 1 cutoff. In this dataset the DPN effects are small: the
#     median |log2FC| of the significant genes is about 0.24 and the largest
#     is about 1.1. A 2-fold cutoff would throw away almost everything.
#     Change lfc_cut below if you want one.
#  b) Three factors, not one. PCA/heatmaps are coloured by patient, time AND
#     treatment, and we add a "cleaned" PCA where patient and time have been
#     removed so any treatment signal becomes visible.
#  c) Two contrasts. DPN is the positive result; OHT is the honest negative
#     one, so we run the same plots on both and compare.
###############################################################################


# =============================================================================
# 0. SETUP
# =============================================================================
# Run these once if a package is missing:
# install.packages(c("ggrepel", "pheatmap", "RColorBrewer", "patchwork"))
# BiocManager::install(c("DESeq2", "apeglm", "clusterProfiler", "enrichplot",
#                        "org.Hs.eg.db", "AnnotationDbi", "limma"))

options(timeout = 1200)

library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)         # non-overlapping gene labels
library(pheatmap)        # heatmaps
library(RColorBrewer)    # colour palettes
library(clusterProfiler) # GO analysis
library(org.Hs.eg.db)    # human annotation
library(enrichplot)      # GO plots
library(patchwork)       # sticking plots together
library(limma)           # removeBatchEffect(), for visualisation only

out_dir <- "analyses/results_DESeq2_parathyroid/visuals_GO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha   <- 0.05   # FDR threshold
lfc_cut <- 0      # minimum |log2FC| to be called up/down (0 = no cutoff, see (a))

# Small helper: save a ggplot into out_dir
save_plot <- function(p, file, width, height) {
  ggsave(file.path(out_dir, file), plot = p, width = width, height = height, dpi = 300)
}


# =============================================================================
# 1. READ THE DATA AND FIT THE MODEL (same as Day 2, model 1b)
# =============================================================================
counts  <- read.csv("Data/parathyroid_counts.csv",  row.names = 1, check.names = FALSE)
samples <- read.csv("Data/parathyroid_samples.csv", row.names = 1)
counts  <- as.matrix(counts)

# patient must be a factor; Control and 24h are the references
samples$patient   <- factor(samples$patient)
samples$treatment <- factor(samples$treatment, levels = c("Control", "DPN", "OHT"))
samples$time      <- factor(samples$time,      levels = c("24h", "48h"))

counts <- counts[, rownames(samples)]
stopifnot(all(colnames(counts) == rownames(samples)))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = samples,
                              design    = ~ patient + time + treatment)

# Same filter as Day 2: >= 10 counts in at least as many samples as the
# smallest treatment group (Control, n = 7)
min_n <- min(table(samples$treatment))
dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
nrow(dds)                       # 16202 genes kept

dds <- DESeq(dds)
resultsNames(dds)

# Gene symbols, looked up once and reused everywhere
symbols <- mapIds(org.Hs.eg.db, keys = rownames(dds), keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
mean(!is.na(symbols))           # fraction of genes that got a name


# =============================================================================
# 2. GENE IDENTIFIER TRANSLATION (needed by clusterProfiler)
#
# clusterProfiler's GO functions want ENTREZ IDs; we have ENSEMBL IDs.
# bitr() = "biological ID translator".
# =============================================================================
id_map <- bitr(rownames(dds), fromType = "ENSEMBL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
head(id_map)

# Not every gene translates. If you lose more than about a third, something
# is wrong with the ID type.
nrow(id_map) / nrow(dds)

# We add ENTREZ IDs with match(), NOT merge(): merge() silently deletes every
# gene without an ENTREZ ID and your counts afterwards come out too low.


# =============================================================================
# 3. RESULT TABLES FOR BOTH CONTRASTS
#
# One helper builds the same tidy data frame for DPN and for OHT:
#   log2FoldChange = SHRUNKEN fold change (for plotting and ranking)
#   lfc_raw        = unshrunk fold change
#   padj           = adjusted p-value from results() (apeglm shrinkage does not
#                    change which genes are significant)
# =============================================================================
make_res_df <- function(dds, trt) {
  res_raw <- results(dds, contrast = c("treatment", trt, "Control"), alpha = alpha)
  res_shr <- lfcShrink(dds, coef = paste0("treatment_", trt, "_vs_Control"),
                       type = "apeglm")
  stopifnot(identical(rownames(res_raw), rownames(res_shr)))

  df <- data.frame(gene           = rownames(res_raw),
                   baseMean       = res_raw$baseMean,
                   log2FoldChange = res_shr$log2FoldChange,   # shrunken
                   lfc_raw        = res_raw$log2FoldChange,
                   stat           = res_raw$stat,
                   pvalue         = res_raw$pvalue,
                   padj           = res_raw$padj,
                   stringsAsFactors = FALSE)

  df$symbol   <- unname(symbols[df$gene])
  df$label    <- ifelse(is.na(df$symbol), df$gene, df$symbol)
  df$ENTREZID <- id_map$ENTREZID[match(df$gene, id_map$ENSEMBL)]

  # Genes with no p-value at all (all-zero counts / outliers) cannot be used
  df <- df[!is.na(df$pvalue), ]

  # Colour category. NA padj (removed by DESeq2's independent filtering)
  # counts as "not significant".
  sig <- !is.na(df$padj) & df$padj < alpha
  df$direction <- "not significant"
  df$direction[sig & df$log2FoldChange >  lfc_cut] <- "up"
  df$direction[sig & df$log2FoldChange < -lfc_cut] <- "down"

  df[order(df$padj), ]
}

res_dpn <- make_res_df(dds, "DPN")
res_oht <- make_res_df(dds, "OHT")

table(res_dpn$direction)   # Day 2 found 135 significant: expect 135 up + down
table(res_oht$direction)   # Day 2 found 13


# =============================================================================
# 4. VARIANCE-STABILISING TRANSFORMATION
#
# Every plot that measures DISTANCE between samples (PCA, heatmaps,
# clustering) needs counts whose variance no longer depends on the mean.
# =============================================================================
vsd     <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
meta    <- as.data.frame(colData(vsd)[, c("patient", "treatment", "time")])


# =============================================================================
# 5. PCA - raw, and with patient + time removed
#
# Raw PCA: what dominates the data. We expect PATIENT and TIME, not drug
# (Day 2 found 8 clusters, one per patient x time).
#
# Cleaned PCA: removeBatchEffect() subtracts the patient and time effects
# while protecting the treatment effect, so what is left is treatment + noise.
# THIS IS FOR PLOTTING ONLY - never feed the cleaned matrix back into DESeq2;
# the model already handles patient and time itself.
# =============================================================================
pca_data    <- plotPCA(vsd, intgroup = c("patient", "treatment", "time"),
                       returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca_raw <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = patient, shape = time)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("Raw VST: patient and time dominate") +
  theme_bw()

design_keep <- model.matrix(~ treatment, data = meta)
vst_clean   <- removeBatchEffect(vst_mat, batch = meta$patient, batch2 = meta$time,
                                 design = design_keep)

pca_clean <- prcomp(t(vst_clean))
pv_clean  <- round(100 * summary(pca_clean)$importance[2, 1:2])
clean_df  <- meta
clean_df$PC1 <- pca_clean$x[, 1]
clean_df$PC2 <- pca_clean$x[, 2]

p_pca_clean <- ggplot(clean_df, aes(x = PC1, y = PC2, colour = treatment, shape = time)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pv_clean[1], "% variance")) +
  ylab(paste0("PC2: ", pv_clean[2], "% variance")) +
  ggtitle("Patient + time removed: treatment visible?") +
  theme_bw()

p_pca <- p_pca_raw + p_pca_clean
p_pca
save_plot(p_pca, "01_PCA_raw_vs_cleaned.png", width = 12, height = 4.8)

# Read it: in the raw plot, colours by patient should form separate groups.
# In the cleaned plot, look for DPN (and OHT) separating from Control. If DPN
# separates but OHT does not, that matches the DE results.


# =============================================================================
# 6. SAMPLE-DISTANCE HEATMAP
#
# Same information as the PCA, as all pairwise distances.
# =============================================================================
sample_dist <- dist(t(vst_mat))        # t(): dist() works on rows, we want samples
dist_matrix <- as.matrix(sample_dist)
rownames(dist_matrix) <- paste(meta$treatment, meta$patient, meta$time, sep = " - ")
colnames(dist_matrix) <- NULL

colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# On screen
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colours, main = "Sample-to-sample distances")
# To file (pheatmap writes straight to the file when filename= is given)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colours, main = "Sample-to-sample distances",
         filename = file.path(out_dir, "02_sample_distances.png"),
         width = 7, height = 5.5)


# =============================================================================
# 7. VOLCANO PLOTS (DPN and OHT, same axes so they can be compared)
# =============================================================================
make_volcano <- function(df, title) {
  d <- df[!is.na(df$padj), ]
  # Label only the 15 most significant genes that pass the threshold
  to_label <- head(d[d$direction != "not significant", ], 15)

  ggplot(d, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = direction), size = 1, alpha = 0.6) +
    scale_colour_manual(values = c("up" = "firebrick",
                                   "down" = "steelblue",
                                   "not significant" = "grey80")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "grey40") +
    geom_text_repel(data = to_label, aes(label = label),
                    size = 3, fontface = "italic", max.overlaps = 20) +
    labs(x = "log2 fold change (shrunken)", y = "-log10 adjusted p-value",
         title = title) +
    theme_bw()
}

p_volc_dpn <- make_volcano(res_dpn, "DPN vs Control")
p_volc_oht <- make_volcano(res_oht, "OHT vs Control")

p_volc_dpn
p_volc_oht
save_plot(p_volc_dpn, "03_volcano_DPN.png", width = 7, height = 5.5)
save_plot(p_volc_oht, "03_volcano_OHT.png", width = 7, height = 5.5)


# =============================================================================
# 8. MA PLOTS (fold change against average expression)
#
# The y-axis range is small here (about +/-1.5), which is the point: these
# are modest shifts, not on/off switches.
# =============================================================================
make_ma <- function(df, title) {
  df$is_sig <- df$direction != "not significant"
  ggplot(df, aes(x = baseMean, y = log2FoldChange, colour = is_sig)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_x_log10() +
    scale_colour_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"),
                        labels = c("TRUE" = "significant", "FALSE" = "not significant"),
                        name = NULL) +
    geom_hline(yintercept = 0) +
    labs(x = "mean of normalised counts", y = "log2 fold change (shrunken)",
         title = title) +
    theme_bw()
}

p_ma <- make_ma(res_dpn, "DPN vs Control") + make_ma(res_oht, "OHT vs Control") +
  plot_layout(guides = "collect")
p_ma
save_plot(p_ma, "04_MA_DPN_OHT.png", width = 12, height = 4.8)


# =============================================================================
# 9. TWO HEATMAPS, AND WHY ONLY ONE OF THEM IS EVIDENCE
# =============================================================================
show_and_save_heatmap <- function(mat, file, title, width = 8, height = 8) {
  pheatmap(mat, annotation_col = meta, main = title)                 # on screen
  pheatmap(mat, annotation_col = meta, main = title,                 # to file
           filename = file.path(out_dir, file), width = width, height = height)
}

# ---- 9a. Top 30 most VARIABLE genes ----------------------------------------
# Chosen without looking at treatment. Expect the samples to split by patient
# and time (that is the real structure in this data), NOT by drug.
gene_var     <- apply(vst_mat, 1, var)
top_variable <- head(order(gene_var, decreasing = TRUE), 30)

mat_var <- vst_mat[top_variable, ]
mat_var <- mat_var - rowMeans(mat_var)                 # center each gene
rownames(mat_var) <- make.unique(ifelse(is.na(symbols[rownames(mat_var)]),
                                        rownames(mat_var),
                                        symbols[rownames(mat_var)]))

show_and_save_heatmap(mat_var, "05_heatmap_variable.png",
                      "Top 30 most variable genes (treatment NOT used to choose them)")

# ---- 9b. Top 30 most SIGNIFICANT DPN genes ---------------------------------
# Chosen BECAUSE they differ with DPN, so this will always look clean:
# circular, for illustration only, never proof that the effect exists.
# We plot the CLEANED matrix so patient and time offsets do not swamp it.
top_sig <- head(res_dpn$gene, 30)

mat_sig <- vst_clean[top_sig, ]
mat_sig <- mat_sig - rowMeans(mat_sig)
rownames(mat_sig) <- make.unique(res_dpn$label[match(top_sig, res_dpn$gene)])

show_and_save_heatmap(mat_sig, "06_heatmap_DPN_significant.png",
                      "Top 30 DPN genes (circular - illustration only; patient/time removed)")


# =============================================================================
# 10. SINGLE-GENE PLOTS (fixes the unlabelled plotCounts from Day 2)
#
# One panel per time point, one line per patient, so you can see whether the
# DPN effect is consistent across all four patients - which is what the
# model assumes.
# =============================================================================
plot_gene <- function(ensg) {
  d <- plotCounts(dds, gene = ensg, intgroup = c("treatment", "time", "patient"),
                  returnData = TRUE, normalized = TRUE)
  sym <- unname(symbols[ensg]); if (is.na(sym)) sym <- ensg

  ggplot(d, aes(x = treatment, y = count, colour = patient, group = patient)) +
    geom_line(alpha = 0.5) +
    geom_point(size = 3) +
    facet_wrap(~ time) +
    scale_y_log10() +
    labs(title = sym, y = "normalised count (log scale)", x = NULL) +
    theme_bw()
}

top_genes <- head(res_dpn$gene, 4)       # the four most significant DPN genes
p_genes <- wrap_plots(lapply(top_genes, plot_gene), ncol = 2) +
  plot_layout(guides = "collect")
p_genes
save_plot(p_genes, "07_top_genes_DPN.png", width = 10, height = 7)


# =============================================================================
# 11. GO OVER-REPRESENTATION ANALYSIS (DPN)
#
# The question: among my significant genes, is any GO term more common than
# you would expect by chance?
#
# THE UNIVERSE MATTERS. The comparison set must be the genes that COULD have
# been detected here, not every gene in the genome. We use every gene that got
# a p-value (i.e. passed the count filter and was tested).
#
# With only 67 up and 68 down genes, the separate lists are small and may give
# few or no terms. So we also run the COMBINED list (all 135), which has the
# most power.
# =============================================================================
universe <- unique(na.omit(res_dpn$ENTREZID))

sig_up   <- res_dpn$ENTREZID[res_dpn$direction == "up"]
sig_down <- res_dpn$ENTREZID[res_dpn$direction == "down"]
sig_all  <- res_dpn$ENTREZID[res_dpn$direction != "not significant"]

length(sig_up); length(sig_down); length(sig_all); length(universe)

# Wrapper that returns NULL (instead of crashing) when nothing is enriched
run_ora <- function(genes, tag) {
  genes <- unique(na.omit(genes))
  if (length(genes) < 5) {
    message(tag, ": fewer than 5 genes, skipped"); return(NULL)
  }
  ego <- enrichGO(gene = genes, universe = universe,
                  OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(tag, ": no enriched GO terms"); return(NULL)
  }
  # GO terms overlap heavily (parents and children); simplify() collapses
  # near-duplicates so the plot is readable.
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

ego_up   <- run_ora(sig_up,   "DPN up")
ego_down <- run_ora(sig_down, "DPN down")
ego_all  <- run_ora(sig_all,  "DPN all")

# How many terms did each list give?
sapply(list(up = ego_up, down = ego_down, all = ego_all),
       function(x) if (is.null(x)) 0 else nrow(as.data.frame(x)))

if (!is.null(ego_all)) {
  head(as.data.frame(ego_all)[, c("Description", "GeneRatio", "p.adjust")], 10)
}

# ---- Dotplots (label_format = 45 keeps long GO names to two lines) ---------
if (!is.null(ego_all)) {
  p_dot_all <- dotplot(ego_all, showCategory = 15, label_format = 45) +
    ggtitle("GO Biological Process: all DPN genes (up + down)")
  print(p_dot_all)
  save_plot(p_dot_all, "08_GO_dotplot_all.png", width = 8.5, height = 7)
}
if (!is.null(ego_up)) {
  p_dot_up <- dotplot(ego_up, showCategory = 12, label_format = 45) +
    ggtitle("GO Biological Process: UP with DPN")
  print(p_dot_up)
  save_plot(p_dot_up, "08_GO_dotplot_up.png", width = 8.5, height = 7)
}
if (!is.null(ego_down)) {
  p_bar_down <- barplot(ego_down, showCategory = 12, label_format = 45) +
    ggtitle("GO Biological Process: DOWN with DPN")
  print(p_bar_down)
  save_plot(p_bar_down, "08_GO_barplot_down.png", width = 8, height = 6.5)
}

# ---- Enrichment map: which terms share genes? ------------------------------
# pairwise_termsim() must be run first. emapplot needs several terms to draw
# a network, so it is wrapped in tryCatch.
if (!is.null(ego_all)) {
  ego_all_sim <- pairwise_termsim(ego_all)
  p_emap <- tryCatch(
    emapplot(ego_all_sim, showCategory = 20) +
      ggtitle("GO terms that share genes (all DPN genes)"),
    error = function(e) { message("emapplot skipped: ", conditionMessage(e)); NULL })
  if (!is.null(p_emap)) {
    print(p_emap)
    save_plot(p_emap, "09_GO_emapplot.png", width = 8, height = 7)
  }
}

# ---- Gene-concept network: which genes drive which term? -------------------
# The plot that turns a term list back into biology. Gene nodes are coloured
# by fold change. readable = TRUE means the object holds SYMBOLS, so the fold
# change vector must be named with symbols too.
if (!is.null(ego_all)) {
  fc_by_symbol <- res_dpn$log2FoldChange
  names(fc_by_symbol) <- res_dpn$symbol
  fc_by_symbol <- fc_by_symbol[!is.na(names(fc_by_symbol))]
  fc_by_symbol <- fc_by_symbol[!duplicated(names(fc_by_symbol))]

  p_cnet <- cnetplot(ego_all, showCategory = 5, foldChange = fc_by_symbol)
  print(p_cnet)
  save_plot(p_cnet, "10_GO_cnetplot.png", width = 9, height = 8)
}


# =============================================================================
# 12. GENE SET ENRICHMENT ANALYSIS (GSEA) - DPN AND OHT
#
# Over-representation needs a cutoff, which throws away everything just below
# it. GSEA uses the WHOLE ranked list, so a pathway where every gene shifts a
# little can still be detected. That makes it the fair test for OHT: is there
# any weak but coordinated response that the 13-gene list cannot show?
#
# We rank by shrunken log2 fold change (as in the template). To rank by the
# Wald statistic instead (accounts for noise per gene), set rank_col <- "stat".
# Here all genes with a fold change are ranked, including those DESeq2's
# independent filtering gave no padj.
# =============================================================================
rank_col <- "log2FoldChange"

run_gsea <- function(df) {
  gl <- df[[rank_col]]
  names(gl) <- df$ENTREZID
  gl <- gl[!is.na(names(gl)) & !is.na(gl)]
  gl <- gl[!duplicated(names(gl))]
  gl <- sort(gl, decreasing = TRUE)          # gseGO requires decreasing order

  set.seed(42)    # GSEA uses permutations: fix the seed for reproducibility
  gseGO(geneList     = gl,
        OrgDb        = org.Hs.eg.db,
        ont          = "BP",
        keyType      = "ENTREZID",
        minGSSize    = 15,
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
}

gsea_dpn <- run_gsea(res_dpn)
gsea_oht <- run_gsea(res_oht)

gsea_dpn_df <- as.data.frame(gsea_dpn)
gsea_oht_df <- as.data.frame(gsea_oht)
nrow(gsea_dpn_df)
nrow(gsea_oht_df)

head(gsea_dpn_df[, c("Description", "setSize", "NES", "p.adjust")], 10)
head(gsea_oht_df[, c("Description", "setSize", "NES", "p.adjust")], 10)

# NES = normalised enrichment score: positive = the set sits at the top of the
# ranking (up with the drug), negative = at the bottom (down).

# Plots, only if something was found
plot_gsea <- function(gsea, tag) {
  if (nrow(as.data.frame(gsea)) == 0) {
    message(tag, ": no enriched gene sets at padj < 0.05"); return(invisible(NULL))
  }
  p_dot <- dotplot(gsea, showCategory = 10, split = ".sign", label_format = 40) +
    facet_grid(. ~ .sign) +
    ggtitle(paste0("GSEA: GO Biological Process (", tag, ")"))
  print(p_dot)
  save_plot(p_dot, paste0("11_GSEA_dotplot_", tag, ".png"), width = 9, height = 7)

  # Ridgeplot: distribution of fold changes inside each set
  p_ridge <- ridgeplot(gsea, showCategory = 15, label_format = 45) +
    ggtitle(paste0("Fold-change distribution within each gene set (", tag, ")"))
  print(p_ridge)
  save_plot(p_ridge, paste0("11_GSEA_ridgeplot_", tag, ".png"), width = 8, height = 7)

  # Running-score plot for the top 3 sets
  p_run <- gseaplot2(gsea, geneSetID = 1:min(3, nrow(as.data.frame(gsea))))
  print(p_run)
  ggsave(file.path(out_dir, paste0("11_GSEA_runningscore_", tag, ".png")),
         plot = p_run, width = 8, height = 6, dpi = 300)
  invisible(NULL)
}

plot_gsea(gsea_dpn, "DPN")
plot_gsea(gsea_oht, "OHT")   # if this is empty, that supports the negative result


# =============================================================================
# 13. COMPOSITE FIGURE (cleaned PCA + DPN volcano + GO dotplot)
# =============================================================================
if (!is.null(ego_all)) {
  p_comp <- (p_pca_clean + p_volc_dpn) / p_dot_all +
    plot_annotation(tag_levels = "A")
  p_comp
  save_plot(p_comp, "12_composite.png", width = 13, height = 11)
}


# =============================================================================
# 14. SAVE THE TABLES
# =============================================================================
write.csv(res_dpn, file.path(out_dir, "DE_results_DPN.csv"), row.names = FALSE)
write.csv(res_oht, file.path(out_dir, "DE_results_OHT.csv"), row.names = FALSE)
if (!is.null(ego_up))   write.csv(as.data.frame(ego_up),   file.path(out_dir, "GO_DPN_up.csv"),   row.names = FALSE)
if (!is.null(ego_down)) write.csv(as.data.frame(ego_down), file.path(out_dir, "GO_DPN_down.csv"), row.names = FALSE)
if (!is.null(ego_all))  write.csv(as.data.frame(ego_all),  file.path(out_dir, "GO_DPN_all.csv"),  row.names = FALSE)
write.csv(gsea_dpn_df, file.path(out_dir, "GSEA_DPN.csv"), row.names = FALSE)
write.csv(gsea_oht_df, file.path(out_dir, "GSEA_OHT.csv"), row.names = FALSE)


# =============================================================================
# 15. THE RESULT IN WORDS
# =============================================================================
n_terms <- function(x) if (is.null(x)) 0 else nrow(as.data.frame(x))

c(DPN_significant        = sum(res_dpn$direction != "not significant"),
  DPN_up                 = sum(res_dpn$direction == "up"),
  DPN_down               = sum(res_dpn$direction == "down"),
  OHT_significant        = sum(res_oht$direction != "not significant"),
  DPN_GO_terms_up        = n_terms(ego_up),
  DPN_GO_terms_down      = n_terms(ego_down),
  DPN_GO_terms_all       = n_terms(ego_all),
  DPN_GSEA_sets          = nrow(gsea_dpn_df),
  OHT_GSEA_sets          = nrow(gsea_oht_df))

sessionInfo()
