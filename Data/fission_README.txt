Dataset: fission
Label  : fission: osmotic stress time course, wt vs atf21 mutant

counts : fission_counts.csv   (7039 genes x 36 samples, raw integers)
samples: fission_samples.csv   (36 rows, 4 columns)

design formula : ~strain + minute
contrast       : minute | 180 | 0
organism       : yeast
gene IDs       : POMBASE
annotation pkg : none for this organism

Load from these files instead of the R package:
    source("scripts/workshop_data.R")
    ds <- load_bulk_dataset("fission", from_csv = TRUE)

Description and your group's task:
Schizosaccharomyces pombe wild type and atf21-deletion mutant sampled 0,
15, 30, 60, 120 and 180 minutes after osmotic stress (1 M sorbitol), in
triplicate (six replicates at time 0). 36 samples - the largest and
best-powered design in the room. You have THREE analyses and the contrast
between them is the lesson: (1) the stress response: 180 min vs 0 min, a
Wald test on the 'minute' coefficient. This is your default contrast and
gives you thousands of genes. (2) the strain main effect, mut vs wt pooled
over time (use contrast_strain). This gives a handful of genes. (3) the
strain:minute interaction, by likelihood-ratio test - the genes whose
RESPONSE TO STRESS differs between strains. This gives fewer than ten
genes. All three numbers are correct. The biology is that osmotic stress
remodels the transcriptome enormously, and that atf21's contribution to
that remodelling is real but very small. Being able to report both without
treating the second as a failure is the skill on display here. Because this
organism has no org.db in our setup, Day 3 enrichment uses a different
route (see the Day 3 script). Reference: Leong et al. 2014 Nat Commun
5:3947.
