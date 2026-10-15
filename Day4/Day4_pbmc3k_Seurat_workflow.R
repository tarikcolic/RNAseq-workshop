###############################################################################
# A COMMON SEURAT WORKFLOW, START TO FINISH
#
# Dataset: pbmc3k - 2,700 peripheral blood mononuclear cells from a healthy
#          donor, 10x Genomics v1 chemistry. Already on the drive, no
#          download needed.
#
# The nine steps every single-cell paper goes through:
#   1  load the matrix          6  cluster the cells
#   2  quality control          7  UMAP
#   3  normalize                8  find marker genes
#   4  find variable genes      9  name the cell types
#   5  scale and run PCA
#
# Figures written to analyses/results_pbmc3k/
#
# Run it line by line with Ctrl+Enter (Windows) or Cmd+Enter (Mac).
# Open 24_Sarajevo_Summer_School.Rproj first, so the file paths below work.
###############################################################################


# =============================================================================
# 0. SETUP
# =============================================================================

library(Seurat)
library(ggplot2)




# =============================================================================
# 1. LOAD THE MATRIX
#
# 10x output is three files in one folder: the counts (matrix.mtx), the gene
# names (genes.tsv) and the cell barcodes (barcodes.tsv). Read10X() reads
# all three and gives you one sparse matrix.
# =============================================================================

counts <- Read10X("/Users/taric/Documents/GitHub/RNAseq-workshop/Data/single_cell/pbmc3k/filtered_gene_bc_matrices/hg19")

# genes x cells
dim(counts)

# Look at a corner of it. The dots are zeros - most of the matrix is zero,
# which is why it is stored as a SPARSE matrix.
counts[1:8, 1:4]

# How much of it is zero?
sum(counts == 0) / length(counts) # 97.41%

# min.cells = 3   keep genes seen in at least 3 cells
# min.features = 200  keep cells with at least 200 genes
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc

# What Seurat already calculated for you:
#   nCount_RNA    total UMIs in that cell
#   nFeature_RNA  how many genes were detected in that cell
head(pbmc@meta.data)



# =============================================================================
# 2. QUALITY CONTROL
#
# Three things make a barcode suspicious:
#   very FEW genes      an empty droplet, or a broken cell
#   very MANY genes     two cells in one droplet (a doublet)
#   high mitochondrial  a dying cell: the cytoplasmic RNA leaked out, so
#                       what is left is disproportionately mitochondrial
# =============================================================================

# Human mitochondrial genes all start with "MT-". In MOUSE they are "mt-".
# Get this wrong and the pattern matches nothing, silently.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Check the pattern actually matched something - 13 genes in human
length(grep("^MT-", rownames(pbmc), value = TRUE))

summary(pbmc$percent.mt)

# Look at the three metrics before choosing any threshold
p_qc <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.1)
p_qc
ggsave("analyses/results_pbmc3k/01_QC_violin.png", plot = p_qc,
       width = 10, height = 4.5, dpi = 300)

# The two relationships worth seeing. Left: cells with high counts but few
# genes are usually debris. Right: the straight line is the expected
# count-gene relationship; points falling off it are the problem cells.
p_sc1 <- FeatureScatter(pbmc, "nCount_RNA", "percent.mt")
p_sc2 <- FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA")
p_scatter <- p_sc1 + p_sc2
p_scatter
ggsave("analyses/results_pbmc3k/02_QC_scatter.png", plot = p_scatter,
       width = 10, height = 4.5, dpi = 300)


# ---- Choose the thresholds -------------------------------------------------
# These three numbers come from the plots above, not from a tutorial.
# For THIS dataset the mitochondrial distribution has a median near 2% and
# a clear tail, so 5% is a reasonable line.
#
# WARNING: 5% is NOT a universal rule. It works here because pbmc3k is old
# v1 chemistry. The same cutoff removes over 97% of the cells in the v3
# datasets in this workshop. Always look at YOUR distribution first.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
                              nFeature_RNA < 2500 &
                              percent.mt   < 5)

# How many cells survived?
pbmc



# =============================================================================
# 3. NORMALISE
#
# A cell with 10,000 UMIs and a cell with 2,000 UMIs are not comparable
# until you divide by the total. LogNormalize does: count / total * 10000,
# then log(x + 1).
# =============================================================================

pbmc <- NormalizeData(pbmc)

# The raw counts are still there - normalization adds a layer, it does not
# overwrite. Compare the same gene before and after:
GetAssayData(pbmc, layer = "counts")[c("CD3D", "LYZ"), 1:5]
GetAssayData(pbmc, layer = "data")[c("CD3D", "LYZ"), 1:5]



# =============================================================================
# 4. FIND VARIABLE GENES
#
# Most genes are the same in every cell and only add noise. Clustering uses
# the genes that actually differ. 2000 is the usual choice.
# =============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# The 10 most variable genes - these should already look like cell type
# markers, which is a good sign
top10 <- head(VariableFeatures(pbmc), 10)
top10

p_var <- LabelPoints(plot = VariableFeaturePlot(pbmc), points = top10,
                     repel = TRUE)
p_var
ggsave("analyses/results_pbmc3k/03_variable_features.png", plot = p_var,
       width = 7.5, height = 5, dpi = 300)



# =============================================================================
# 5. SCALE AND RUN PCA
#
# ScaleData centres each gene on zero and sets its variance to 1, so that a
# highly expressed gene does not dominate the PCA just by being loud.
# =============================================================================

pbmc <- ScaleData(pbmc)

# PCA runs on the variable genes by default
pbmc <- RunPCA(pbmc)

# What drives the first PCs?
print(pbmc[["pca"]], dims = 1:3, nfeatures = 8)

DimPlot(pbmc, reduction = "pca")

# DimHeatmap shows the cells and genes at both extremes of a PC. A PC that
# looks like clean structure is real; one that looks like noise is noise.
DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()


# ---- How many PCs should you keep? -----------------------------------------
# The elbow plot: find where the curve flattens. Everything after that is
# mostly noise. There is no exact answer - 10 is a defensible read here.
p_elbow <- ElbowPlot(pbmc, ndims = 40)
p_elbow
ggsave("analyses/results_pbmc3k/06_elbow.png", plot = p_elbow,
       width = 6, height = 4.5, dpi = 300)

# Being slightly too generous is safer than being too strict.
n_pcs <- 12



# =============================================================================
# 6. CLUSTER THE CELLS
#
# Two steps. FindNeighbors builds a graph joining each cell to the cells
# most like it. FindClusters cuts that graph into communities.
#
# resolution controls how many clusters you get. Higher = more clusters.
# 0.4 to 1.2 is the usual range for a dataset this size.
# =============================================================================

pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# How many clusters, and how big?
table(Idents(pbmc))

# The resolution is a CHOICE, not a result. Try a few and see what moves.



# Go back to the one we will use
Idents(pbmc) <- pbmc$seurat_clusters
table(Idents(pbmc))



# =============================================================================
# 7. UMAP
#
# UMAP is for LOOKING at the data. Distances between clusters on a UMAP are
# not meaningful, and neither is their arrangement. The clustering was done
# on the PCs, not on the UMAP.
# =============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)

p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5) +
  ggtitle("pbmc3k")
p_umap
ggsave("analyses/results_pbmc3k/07_UMAP_clusters.png", plot = p_umap,
       width = 7, height = 5.5, dpi = 300)

# Always check the QC metrics on the UMAP too. A cluster that is just
# "cells with high mitochondrial content" is a technical artefact, not a
# cell type.
p_umap_qc <- FeaturePlot(pbmc, features = c("nFeature_RNA", "percent.mt"))
p_umap_qc
ggsave("analyses/results_pbmc3k/08_UMAP_QC.png", plot = p_umap_qc,
       width = 10, height = 4.5, dpi = 300)



# =============================================================================
# 8. FIND MARKER GENES
#
# For each cluster: which genes are higher in it than in all other cells?
#
# only.pos = TRUE      we want markers OF the cluster, not against it
# min.pct = 0.25       the gene must be detected in 25% of cells
# logfc.threshold      ignore tiny differences, for speed
# =============================================================================

markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)

dim(markers)
head(markers)

# The top 3 markers per cluster, by fold change
markers_sorted <- markers[order(markers$cluster, -markers$avg_log2FC), ]
top3 <- do.call(rbind, by(markers_sorted, markers_sorted$cluster, head, 3))
top3[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]

write.csv(markers, "analyses/results_pbmc3k/cluster_markers.csv",
          row.names = FALSE)


# ---- Look at the markers three ways ----------------------------------------
# On the UMAP
p_feat <- FeaturePlot(pbmc, features = c("MS4A1", "CD3D", "LYZ", "NKG7"))
p_feat
ggsave("analyses/results_pbmc3k/09_markers_umap.png", plot = p_feat,
       width = 9, height = 7, dpi = 300)

# As distributions per cluster - this shows you whether a marker is
# genuinely specific or just slightly higher
p_vln <- VlnPlot(pbmc, features = c("MS4A1", "CD3D", "LYZ"), ncol = 3)
p_vln
ggsave("analyses/results_pbmc3k/10_markers_violin.png", plot = p_vln,
       width = 11, height = 4, dpi = 300)

# As a heatmap of the top 10 per cluster
markers_top10 <- do.call(rbind, by(markers_sorted, markers_sorted$cluster,
                                   head, 10))
p_heat <- DoHeatmap(pbmc, features = markers_top10$gene) + NoLegend()
p_heat
ggsave("analyses/results_pbmc3k/11_marker_heatmap.png", plot = p_heat,
       width = 10, height = 11, dpi = 300)



# =============================================================================
# 9. NAME THE CELL TYPES
#
# This is the step that needs biology, not code. You match each cluster's
# markers against genes you know.
#
# These are the canonical PBMC markers:
#   IL7R, CCR7     naive CD4 T
#   IL7R, S100A4   memory CD4 T
#   CD8A           CD8 T
#   MS4A1          B
#   CD14, LYZ      CD14+ monocyte
#   FCGR3A, MS4A7  FCGR3A+ monocyte
#   GNLY, NKG7     NK
#   FCER1A, CST3   dendritic cell
#   PPBP           platelet
# =============================================================================

# A dotplot is the fastest way to see which cluster has which marker.
# Dot size = fraction of cells expressing, colour = average level.
canonical <- c("IL7R", "CCR7", "S100A4", "CD8A", "MS4A1", "CD14", "LYZ",
               "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")

p_dot <- DotPlot(pbmc, features = canonical) + RotatedAxis()
p_dot
ggsave("analyses/results_pbmc3k/12_canonical_dotplot.png", plot = p_dot,
       width = 9, height = 5, dpi = 300)

# Read the dotplot, then write the names down IN CLUSTER ORDER (0, 1, 2...).
# Check this against YOUR dotplot - cluster numbering is not guaranteed to
# be stable if you change the resolution or the number of PCs.
new_names <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B",
               "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

# There must be one name per cluster
length(new_names) == length(levels(pbmc))

names(new_names) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_names)

p_final <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4,
                   repel = TRUE) + NoLegend() +
  ggtitle("pbmc3k: annotated cell types")
p_final
ggsave("analyses/results_pbmc3k/13_UMAP_annotated.png", plot = p_final,
       width = 7, height = 5.5, dpi = 300)

# How many cells of each type?
table(Idents(pbmc))



# =============================================================================
# 10. SAVE
# =============================================================================

# Store the labels in the metadata so they survive a save and reload
pbmc$cell_type <- Idents(pbmc)

saveRDS(pbmc, "analyses/results_pbmc3k/pbmc3k_annotated.rds")

# To pick up where you left off next time:
#   pbmc <- readRDS("analyses/results_pbmc3k/pbmc3k_annotated.rds")

figure <- (p_umap | p_final) / p_dot
ggsave("analyses/results_pbmc3k/14_composite.png", plot = figure,
       width = 13, height = 10, dpi = 300)


# =============================================================================
# 11. THE RESULT IN WORDS
# =============================================================================

ncol(pbmc)                      # cells kept
nrow(pbmc)                      # genes kept
length(levels(pbmc$cell_type))  # cell types found
table(pbmc$cell_type)

sessionInfo()
