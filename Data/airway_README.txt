Dataset: airway
Label  : airway: dexamethasone vs untreated, 4 donors

counts : airway_counts.csv   (63677 genes x 8 samples, raw integers)
samples: airway_samples.csv   (8 rows, 9 columns)

design formula : ~cell + dex
contrast       : dex | trt | untrt
organism       : human
gene IDs       : ENSEMBL
annotation pkg : org.Hs.eg.db

Load from these files instead of the R package:
    source("scripts/workshop_data.R")
    ds <- load_bulk_dataset("airway", from_csv = TRUE)

Description and your group's task:
Primary human airway smooth muscle cells from four donors were treated with
1 uM dexamethasone for 18 h or left untreated. Because each donor
contributes both a treated and an untreated sample, this is a paired
design: include 'cell' in the model to absorb between-donor variation.
Biological question: which genes does a glucocorticoid regulate in airway
smooth muscle? Reference: Himes et al. 2014 PLoS One 9(6):e99625; GEO
GSE52778.
