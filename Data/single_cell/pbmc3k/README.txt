==============================================================================
pbmc3k  —  GROUP 1
PBMC 3k (human, 10x v1 chemistry, 2,700 cells)
==============================================================================

SOURCE
  10x Genomics — 3k PBMCs from a healthy donor
  https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0
  Public 10x Genomics reference dataset. Redistributed here unmodified;
  these are the files as downloaded from 10x, not reprocessed.

FILES
  pbmc3k/filtered_gene_bc_matrices/
    hg19
  0.0 MB on disk. Standard 10x Market Exchange (MEX) triplet.

LOAD IT
  source("scripts/workshop_sc_data.R")
  sc  <- load_sc_dataset("pbmc3k")       # prints this dataset's task
  obj <- Seurat::CreateSeuratObject(Seurat::Read10X(sc$path),
                                    project = "pbmc3k",
                                    min.cells = 3, min.features = 200)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = sc$mt_pattern)

  The files are already on the drive. You do NOT need to download anything.

WHAT IS IN IT  (measured from the matrix in this folder, before any filtering)
  organism                 human
  chemistry                v1 (2016)
  barcodes (cells)         2,700
  features (genes)         32,738   of which 16,634 are detected at all
  matrix sparsity          97.41% zeros

  UMIs per cell            median 2,197
                           5-95%  956 - 4,221
                           IQR    1,758 - 2,763
  genes per cell           median 817
                           5-95%  436 - 1,368

  mitochondrial content    pattern "^MT-" matches 13 genes
                           median 2.03%   75th 2.64%   95th 4.01%
                           99th 5.88%   max 22.57%
  ribosomal content        pattern "^RP[SL]" matches 106 genes
                           median 36.77%

  HOW MANY CELLS A MITOCHONDRIAL THRESHOLD WOULD REMOVE
                           > 5%  ->  2.1% of cells
                           >10%  ->  0.2% of cells
                           >20%  ->  0% of cells

  Suggested STARTING POINT for this dataset: percent.mt < 5
  That is a starting point, not a rule. You must justify whatever you use.

MOST ABUNDANT GENES (share of all counts in the matrix)
  `MALAT1` (2.53%), `TMSB4X` (1.94%), `B2M` (1.9%), `RPL10` (1.39%), `RPL13` (1.21%), `RPL13A` (1.2%)

CANONICAL MARKERS FOR ANNOTATION (Day 5)
| Expected population | Gene |
|---|---|
| Naive CD4 T | `IL7R` |
| Memory CD4 T | `S100A4` |
| CD8 T | `CD8A` |
| NK | `GNLY` |
| B | `MS4A1` |
| CD14+ Mono | `CD14` |
| FCGR3A+ Mono | `FCGR3A` |
| DC | `FCER1A` |
| Platelet | `PPBP` |

YOUR TASK
  Peripheral blood mononuclear cells from a healthy donor, the reference
  dataset of the single-cell field and the one every tutorial uses. 10x v1
  chemistry, so UMI counts per cell are LOW by modern standards (~2,000
  median) - a useful contrast with group 2, who have the same tissue on v3
  chemistry. Cell types present: CD4 and CD8 T cells, NK cells, B cells,
  CD14+ and FCGR3A+ monocytes, dendritic cells, platelets. Expect 8-9
  clusters. Your task at the round-up: how many UMIs and genes per cell do
  you see, and how does that compare with group 2?
