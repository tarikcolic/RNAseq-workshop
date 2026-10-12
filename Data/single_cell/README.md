# Single-cell datasets

Five public 10x Genomics reference datasets, one per group, already
downloaded. **Nothing here needs an internet connection during the workshop.**

```
datasets/single_cell/
├── README.md                    <- this file
├── pbmc3k/      README.txt + filtered_gene_bc_matrices/hg19/{barcodes,genes,matrix}
├── pbmc1k_v3/   README.txt + filtered_feature_bc_matrix/{barcodes,features,matrix}.gz
├── neuron1k_v3/ README.txt + filtered_feature_bc_matrix/...
├── pbmc5k_v3/   README.txt + filtered_feature_bc_matrix/...
└── heart1k_v3/  README.txt + filtered_feature_bc_matrix/...
```

Each folder has its own `README.txt` with that dataset's full QC profile,
canonical markers, and the group's specific task. Load any of them with:

```r
source("scripts/workshop_sc_data.R")
sc  <- load_sc_dataset("pbmc3k")                 # prints the task
obj <- Seurat::CreateSeuratObject(Seurat::Read10X(sc$path),
                                  min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = sc$mt_pattern)
```

## The five, side by side

All numbers measured from the matrices in this folder, before any filtering.

| Group | Dataset | Organism | Chem | Cells | Median UMI | Median genes | Median %MT | Cells lost at %MT < 5 |
|---|---|---|---|---|---|---|---|---|
| 1 | `pbmc3k` | human | v1 (2016) | 2,700 | 2,197 | 817 | 2.03% | **2.1%** |
| 2 | `pbmc1k_v3` | human | v3 | 1,222 | 6,628 | 1,919 | 10.29% | **98.2%** |
| 3 | `neuron1k_v3` | mouse | v3 | 1,301 | 9,744 | 3,356 | 7.16% | **81.5%** |
| 4 | `pbmc5k_v3` | human | v3 | 5,025 | 7,803 | 2,171 | 8.45% | **97.3%** |
| 5 | `heart1k_v3` | mouse | v3 | 1,011 | 7,610 | 2,058 | 7.21% | **73.5%** |

## Read that last column before you teach this

The Seurat PBMC tutorial filters `percent.mt < 5`. Every student will have
seen it, and most will copy it. Applied to these datasets it removes:

- `pbmc3k` — **2.1%** of cells. Harmless.
- `heart1k_v3` — **73.5%**
- `neuron1k_v3` — **81.5%**
- `pbmc5k_v3` — **97.3%**
- `pbmc1k_v3` — **98.2%**

**The threshold is an artifact of v1 chemistry, not a biological rule.** The
tutorial was written on `pbmc3k`, a v1 dataset with a median of
2.03% mitochondrial content. v3 chemistry captures far more
mitochondrial transcript — median 10.29% in `pbmc1k_v3` —
so the same number applied to a v3 dataset destroys it. The two v3 **PBMC**
datasets are hit hardest of all: same tissue, same species, same healthy
donors as the tutorial, and the threshold still deletes ~97-98% of the cells.

This is the single most useful thing in the single-cell half of the week,
because it generalises: a parameter published with one dataset is not a
constant of nature, and the failure is silent. Nobody's code errors. The UMAP
still looks fine. There are just 22 cells in it.

Where the tissues genuinely differ is in the **tail**, not the median.
`heart1k_v3` has 24.9% of cells above 20%
mitochondrial content, against ~9-12% for the others — that is real
cardiomyocyte biology on top of the chemistry effect, and the reason Group 5
has the hardest threshold call to defend.

## Two other things visible in the numbers

**Dissociated solid tissue carries blood.** `Hbb-bs` is the single most
abundant gene in `heart1k_v3` (4.42% of all
counts) and is in the top five of `neuron1k_v3`. Groups 3 and 5 should expect
an erythrocyte cluster and decide what to do about it. The PBMC datasets show
`MALAT1` on top instead, which is nuclear and normal.

**Depth is not comparable across chemistries.** `pbmc3k` has a median of
2,197 UMIs and 817 genes per
cell; the v3 PBMC sets have roughly 3x the UMIs and 2.5x the genes. When
Group 1 and Group 2 compare their UMAPs at the round-up, they are not looking
at the same resolution of data. That difference is chemistry, not biology.

## Sizes

| Dataset | On disk |
|---|---|
| `pbmc3k` | 27.7 MB |
| `pbmc1k_v3` | 9.2 MB |
| `neuron1k_v3` | 14.9 MB |
| `pbmc5k_v3` | 41.1 MB |
| `heart1k_v3` | 9.5 MB |
| **total** | **102.4 MB** |

`pbmc3k` is the largest despite having the fewest UMIs because its files are
uncompressed — it predates the gzipped MEX convention. Seurat reads both.
