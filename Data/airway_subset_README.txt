Dataset: airway_subset
Label  : airway (3 v 3, unpaired): dexamethasone vs untreated

counts : airway_subset_counts.csv   (63677 genes x 6 samples, raw integers)
samples: airway_subset_samples.csv   (6 rows, 9 columns)

design formula : ~dex
contrast       : dex | trt | untrt
organism       : human
gene IDs       : ENSEMBL
annotation pkg : org.Hs.eg.db

Load from these files instead of the R package:
    source("scripts/workshop_data.R")
    ds <- load_bulk_dataset("airway_subset", from_csv = TRUE)

Description and your group's task:
The airway experiment reduced to three donors and analysed as an UNPAIRED 3
v 3 comparison, ignoring the donor pairing. Your dataset is the control
experiment for the whole room. Group 1 has the same biology with four
donors and a paired design. Expect to find substantially fewer significant
genes than they do. Your task at the round-up is to explain WHY, and to
demonstrate it: re-run with design = ~ cell + dex and compare the number of
hits. This is the clearest demonstration you will ever see of what
experimental design buys you. Reference: Himes et al. 2014 PLoS One
9(6):e99625; GEO GSE52778.
