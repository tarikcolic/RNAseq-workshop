Dataset: parathyroid
Label  : parathyroid: DPN vs Control, 4 patients, 24/48 h

counts : parathyroid_counts.csv   (63193 genes x 23 samples, raw integers)
samples: parathyroid_samples.csv   (23 rows, 3 columns)

design formula : ~patient + time + treatment
contrast       : treatment | DPN | Control
organism       : human
gene IDs       : ENSEMBL
annotation pkg : org.Hs.eg.db

Load from these files instead of the R package:
    source("scripts/workshop_data.R")
    ds <- load_bulk_dataset("parathyroid", from_csv = TRUE)

Description and your group's task:
Primary cultures of human parathyroid adenoma from four patients, treated
with the oestrogen-receptor-beta agonist DPN, with 4-hydroxytamoxifen
(OHT), or with vehicle control, harvested at 24 h and 48 h. Technical
replicates have been summed, leaving 23 libraries. This is the only
genuinely three-factor design in the room (patient x treatment x time), and
you have two tasks: (1) DOES THE MODEL MATTER? Run DPN vs Control with ~
patient + treatment, then with ~ patient + time + treatment. The number of
significant genes changes by roughly eight-fold from the same reads.
Nothing about the biology changed; you simply stopped letting a large
time-in-culture effect masquerade as noise. Report both numbers and explain
which you would publish. (2) A HONEST NEGATIVE RESULT. Now run OHT vs
Control. You will find very little. That is the answer, not a failure - and
writing it up correctly ('we detected no reproducible transcriptional
response to OHT at this power') is a harder and more useful skill than
finding a long gene list. Biological question: does ER-beta signalling
regulate gene expression in parathyroid tumour cells? Reference: Haglund et
al. 2012 JCEM 97:4631; GEO GSE37211.
