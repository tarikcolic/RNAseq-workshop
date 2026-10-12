==============================================================================
neuron1k_v3  —  GROUP 3
Neurons 1k (mouse E18 cortex, 10x v3, 1,301 cells)
==============================================================================

SOURCE
  10x Genomics — 1k brain cells from an E18 mouse (v3)
  https://www.10xgenomics.com/datasets/1-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0
  Public 10x Genomics reference dataset. Redistributed here unmodified;
  these are the files as downloaded from 10x, not reprocessed.

FILES
  neuron1k_v3/filtered_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  14.9 MB on disk. Standard 10x Market Exchange (MEX) triplet.

LOAD IT
  source("scripts/workshop_sc_data.R")
  sc  <- load_sc_dataset("neuron1k_v3")       # prints this dataset's task
  obj <- Seurat::CreateSeuratObject(Seurat::Read10X(sc$path),
                                    project = "neuron1k_v3",
                                    min.cells = 3, min.features = 200)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = sc$mt_pattern)

  The files are already on the drive. You do NOT need to download anything.

WHAT IS IN IT  (measured from the matrix in this folder, before any filtering)
  organism                 mouse
  chemistry                v3
  barcodes (cells)         1,301
  features (genes)         31,053   of which 18,085 are detected at all
  matrix sparsity          89.55% zeros

  UMIs per cell            median 9,744
                           5-95%  1,013 - 25,885
                           IQR    5,333 - 15,821
  genes per cell           median 3,356
                           5-95%  238 - 5,640

  mitochondrial content    pattern "^mt-" matches 13 genes
                           median 7.16%   75th 11.41%   95th 46%
                           99th 88.77%   max 97.53%
  ribosomal content        pattern "^Rp[sl]" matches 105 genes
                           median 11.5%

  HOW MANY CELLS A MITOCHONDRIAL THRESHOLD WOULD REMOVE
                           > 5%  ->  81.5% of cells
                           >10%  ->  29.7% of cells
                           >20%  ->  12.4% of cells

  Suggested STARTING POINT for this dataset: percent.mt < 10
  That is a starting point, not a rule. You must justify whatever you use.

MOST ABUNDANT GENES (share of all counts in the matrix)
  `mt-Atp6` (1.67%), `Malat1` (1.61%), `mt-Co3` (1.51%), `Hbb-bs` (1.18%), `mt-Co2` (1.07%), `Tuba1a` (1.01%)

CANONICAL MARKERS FOR ANNOTATION (Day 5)
| Expected population | Gene |
|---|---|
| Radial glia | `Sox2` |
| Cycling | `Mki67` |
| Intermediate prog | `Eomes` |
| Immature neuron | `Neurod2` |
| Excitatory neuron | `Neurod6` |
| Interneuron | `Gad1` |
| Oligodendrocyte prog | `Olig1` |
| Microglia | `Cx3cr1` |
| Endothelial | `Cldn5` |

YOUR TASK
  Cells dissociated from the cortex of an embryonic day 18 mouse. A
  developing tissue, so you should expect a CONTINUUM - radial glia,
  intermediate progenitors, immature neurons - rather than discrete types.
  This is the dataset in the room where clustering is most obviously a
  decision rather than a discovery. NOTE: mouse mitochondrial genes are
  'mt-', not 'MT-'. Your task at the round-up: show your UMAP at three
  resolutions and argue for one. Is there a correct number of clusters here?
