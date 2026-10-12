Dataset: pasilla
Label  : pasilla: RNAi knockdown vs untreated, D. melanogaster S2

counts : pasilla_counts.csv   (14599 genes x 7 samples, raw integers)
samples: pasilla_samples.csv   (7 rows, 2 columns)

design formula : ~type + condition
contrast       : condition | treated | untreated
organism       : fly
gene IDs       : FLYBASE
annotation pkg : org.Dm.eg.db

Load from these files instead of the R package:
    source("scripts/workshop_data.R")
    ds <- load_bulk_dataset("pasilla", from_csv = TRUE)

Description and your group's task:
Drosophila S2 cells with RNAi knockdown of the splicing factor pasilla (3
treated) versus untreated controls (4 samples). Libraries were prepared as
a mixture of single-end and paired-end runs, so 'type' is a technical
covariate that must be in the model. Biological question: what does loss of
a splicing regulator do to steady-state gene expression? Note that the
primary effect of pasilla loss is on SPLICING, which gene-level counts
capture only indirectly - a good discussion point. Reference: Brooks et al.
2011 Genome Res 21:193; GEO GSE18508.
