==============================================================================
heart1k_v3  —  GROUP 5
Heart 1k (mouse E18 heart, 10x v3, 1,011 cells)
==============================================================================

SOURCE
  10x Genomics — 1k heart cells from an E18 mouse (v3)
  https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0
  Public 10x Genomics reference dataset. Redistributed here unmodified;
  these are the files as downloaded from 10x, not reprocessed.

FILES
  heart1k_v3/filtered_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  9.5 MB on disk. Standard 10x Market Exchange (MEX) triplet.

LOAD IT
  source("scripts/workshop_sc_data.R")
  sc  <- load_sc_dataset("heart1k_v3")       # prints this dataset's task
  obj <- Seurat::CreateSeuratObject(Seurat::Read10X(sc$path),
                                    project = "heart1k_v3",
                                    min.cells = 3, min.features = 200)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = sc$mt_pattern)

  The files are already on the drive. You do NOT need to download anything.

WHAT IS IN IT  (measured from the matrix in this folder, before any filtering)
  organism                 mouse
  chemistry                v3
  barcodes (cells)         1,011
  features (genes)         31,053   of which 19,436 are detected at all
  matrix sparsity          91.56% zeros

  UMIs per cell            median 7,610
                           5-95%  789 - 46,641
                           IQR    2,846 - 19,775
  genes per cell           median 2,058
                           5-95%  64 - 6,357

  mitochondrial content    pattern "^mt-" matches 13 genes
                           median 7.21%   75th 19.96%   95th 50.47%
                           99th 82.98%   max 98.96%
  ribosomal content        pattern "^Rp[sl]" matches 105 genes
                           median 14.05%

  HOW MANY CELLS A MITOCHONDRIAL THRESHOLD WOULD REMOVE
                           > 5%  ->  73.5% of cells
                           >10%  ->  35.7% of cells
                           >20%  ->  24.9% of cells

  Suggested STARTING POINT for this dataset: percent.mt < 30
  That is a starting point, not a rule. You must justify whatever you use.

MOST ABUNDANT GENES (share of all counts in the matrix)
  `Hbb-bs` (4.42%), `mt-Co3` (1.78%), `mt-Atp6` (1.68%), `mt-Co1` (1.61%), `mt-Co2` (1.4%), `Malat1` (1.21%)

CANONICAL MARKERS FOR ANNOTATION (Day 5)
| Expected population | Gene |
|---|---|
| Cardiomyocyte | `Actc1` |
| Cardiomyocyte 2 | `Tnnt2` |
| Fibroblast | `Col1a1` |
| Endothelial | `Cdh5` |
| Endocardial | `Npr3` |
| Epicardial | `Wt1` |
| Smooth muscle | `Acta2` |
| Immune | `Ptprc` |
| Erythrocyte | `Hba-a1` |

YOUR TASK
  Cells from an embryonic day 18 mouse heart, and the dataset with the
  hardest threshold call in the room. The percent.mt < 5 cutoff everyone
  copies from the Seurat PBMC tutorial removes 73.5% of these cells. The
  obvious conclusion is 'cardiomyocytes are mitochondria-rich, so heart
  needs a special threshold'. MEASURE IT BEFORE YOU BELIEVE IT: the same
  cutoff removes 98.2% of pbmc1k_v3 and 97.3% of pbmc5k_v3 - healthy human
  blood, the very tissue the tutorial was written on. The cutoff is an
  artefact of the OLD v1 chemistry used for pbmc3k (median 2.0%
  mitochondrial); v3 captures far more mitochondrial transcript, and every
  v3 dataset in this room is wrecked by it. What IS specific to your heart
  tissue is the TAIL rather than the median: 24.9% of these cells sit above
  20% mitochondrial content, against 9-12% in the other v3 datasets. Some of
  those are dying cells you should remove; some are cardiomyocytes you must
  keep, and no threshold separates them cleanly. That is the real problem,
  and 'the tutorial said 5' is not an answer to it. Your task at the round-
  up: show your percent.mt distribution, state the threshold you chose and
  the reasoning, and show what the UMAP looks like at 5% instead. NOTE:
  mouse genes are 'mt-', not 'MT-'. Two more things the counts already tell
  you: Hbb-bs is the single most abundant gene here (4.4% of all counts), so
  expect an erythrocyte cluster from blood in the dissociated tissue; and as
  a solid tissue this was enzymatically dissociated, so check for the stress
  signature (Fos, Jun, Hspa1a) before believing any cluster defined by it.
