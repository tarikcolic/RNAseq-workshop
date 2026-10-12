==============================================================================
pbmc5k_v3  —  GROUP 4
PBMC 5k (human, 10x v3, 5,025 cells)
==============================================================================

SOURCE
  10x Genomics — 5k peripheral blood mononuclear cells (v3)
  https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-1-0
  Public 10x Genomics reference dataset. Redistributed here unmodified;
  these are the files as downloaded from 10x, not reprocessed.

FILES
  pbmc5k_v3/filtered_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  41.1 MB on disk. Standard 10x Market Exchange (MEX) triplet.

LOAD IT
  source("scripts/workshop_sc_data.R")
  sc  <- load_sc_dataset("pbmc5k_v3")       # prints this dataset's task
  obj <- Seurat::CreateSeuratObject(Seurat::Read10X(sc$path),
                                    project = "pbmc5k_v3",
                                    min.cells = 3, min.features = 200)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = sc$mt_pattern)

  The files are already on the drive. You do NOT need to download anything.

WHAT IS IN IT  (measured from the matrix in this folder, before any filtering)
  organism                 human
  chemistry                v3
  barcodes (cells)         5,025
  features (genes)         33,538   of which 21,819 are detected at all
  matrix sparsity          93.25% zeros

  UMIs per cell            median 7,803
                           5-95%  1,566 - 17,867
                           IQR    6,079 - 10,210
  genes per cell           median 2,171
                           5-95%  597 - 4,008

  mitochondrial content    pattern "^MT-" matches 13 genes
                           median 8.45%   75th 10.54%   95th 35.49%
                           99th 85.7%   max 96.79%
  ribosomal content        pattern "^RP[SL]" matches 104 genes
                           median 27.97%

  HOW MANY CELLS A MITOCHONDRIAL THRESHOLD WOULD REMOVE
                           > 5%  ->  97.3% of cells
                           >10%  ->  29.3% of cells
                           >20%  ->  9.7% of cells

  Suggested STARTING POINT for this dataset: percent.mt < 15
  That is a starting point, not a rule. You must justify whatever you use.

MOST ABUNDANT GENES (share of all counts in the matrix)
  `MALAT1` (4.11%), `MT-CO1` (1.66%), `MT-CO2` (1.24%), `MT-CO3` (1.2%), `EEF1A1` (1.17%), `MT-ATP6` (1.09%)

CANONICAL MARKERS FOR ANNOTATION (Day 5)
| Expected population | Gene |
|---|---|
| CD4 T | `IL7R` |
| CD8 T | `CD8A` |
| NK | `GNLY` |
| B | `MS4A1` |
| CD14+ Mono | `CD14` |
| FCGR3A+ Mono | `FCGR3A` |
| cDC | `FCER1A` |
| pDC | `LILRA4` |
| Platelet | `PPBP` |
| Plasma | `JCHAIN` |

YOUR TASK
  Peripheral blood mononuclear cells, 5,000 cells on v3 chemistry. You have
  the largest dataset in the room, which means two things: your rare
  populations (dendritic cells, plasmacytoid DCs, megakaryocytes) are
  actually resolvable, and your DOUBLET rate is the highest - around 4%,
  i.e. roughly 200 of your 'cells' are two cells. Your task at the round-up:
  find the doublet cluster. Look for a cluster with unusually high
  nCount_RNA that co-expresses two lineage programmes (e.g. CD3D and LYZ, or
  MS4A1 and CD14).
