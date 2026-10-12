==============================================================================
pbmc1k_v3  —  GROUP 2
PBMC 1k (human, 10x v3 chemistry, 1,222 cells)
==============================================================================

SOURCE
  10x Genomics — 1k PBMCs from a healthy donor (v3)
  https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
  Public 10x Genomics reference dataset. Redistributed here unmodified;
  these are the files as downloaded from 10x, not reprocessed.

FILES
  pbmc1k_v3/filtered_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  9.2 MB on disk. Standard 10x Market Exchange (MEX) triplet.

LOAD IT
  source("scripts/workshop_sc_data.R")
  sc  <- load_sc_dataset("pbmc1k_v3")       # prints this dataset's task
  obj <- Seurat::CreateSeuratObject(Seurat::Read10X(sc$path),
                                    project = "pbmc1k_v3",
                                    min.cells = 3, min.features = 200)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = sc$mt_pattern)

  The files are already on the drive. You do NOT need to download anything.

WHAT IS IN IT  (measured from the matrix in this folder, before any filtering)
  organism                 human
  chemistry                v3
  barcodes (cells)         1,222
  features (genes)         33,538   of which 18,391 are detected at all
  matrix sparsity          93.9% zeros

  UMIs per cell            median 6,628
                           5-95%  1,026 - 15,945
                           IQR    4,907 - 8,965
  genes per cell           median 1,919
                           5-95%  448 - 3,740

  mitochondrial content    pattern "^MT-" matches 13 genes
                           median 10.29%   75th 12.59%   95th 46.86%
                           99th 98.51%   max 99.7%
  ribosomal content        pattern "^RP[SL]" matches 104 genes
                           median 25.84%

  HOW MANY CELLS A MITOCHONDRIAL THRESHOLD WOULD REMOVE
                           > 5%  ->  98.2% of cells
                           >10%  ->  54.6% of cells
                           >20%  ->  8.8% of cells

  Suggested STARTING POINT for this dataset: percent.mt < 15
  That is a starting point, not a rule. You must justify whatever you use.

MOST ABUNDANT GENES (share of all counts in the matrix)
  `MALAT1` (3.76%), `MT-CO1` (2.01%), `MT-CO3` (1.86%), `MT-CO2` (1.45%), `MT-ATP6` (1.26%), `EEF1A1` (1.21%)

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
  Peripheral blood mononuclear cells from a healthy donor on 10x v3
  chemistry. Same tissue as group 1, much higher capture efficiency: expect
  ~5,000 median UMIs per cell against their ~2,000. Fewer cells, better data
  per cell. Your task at the round-up: compare your QC distributions with
  group 1's and decide whether the extra sensitivity changed the biology you
  could see, or only the numbers. Do you resolve any cell type they cannot?
