# Ghalambor 2015 – Guppy brain RNA-seq

**Paper:** Ghalambor CK, Hoke KL, Ruell EW, Fischer EK, Reznick DN, Hughes KA (2015). Non-adaptive plasticity potentiates rapid adaptive evolution of gene expression in nature. *Nature* 525:372–375. https://doi.org/10.1038/nature15256

**Data source:** NCBI BioProject PRJNA291937 (SRA study SRP062364)

## The data

- **Species:** Trinidadian guppy (*Poecilia reticulata*)
- **Tissue:** whole brain from adult males. Each library pools the brains of two full siblings from the same treatment. Reads are single-end 100 bp, sequenced on an Illumina HiSeq 2000.
- **Populations** (the prefix on the sample names):

| Code | Population | Role |
|---|---|---|
| `gh` | Guanapo, high-predation | **Ancestral source population**: lives with predatory cichlids |
| `ul` | Upper Lalaja | **Introduced** population: Guanapo fish moved to a cichlid-free stream |
| `ll` | Lower Lalaja | **Introduced** population: a second, independent replicate of the same transplant |
| `tl` | Taylor, low-predation | **Native** cichlid-free population: shows what long-term adaptation to low predation looks like |

- **Rearing treatment** (the `S1`/`S2` part of the sample names): the fish were raised in the lab either **without** predator cues (`S1` = "NP") or **with** predator cues (`S2` = "P"). `R1`–`R5` are the replicate libraries. This gives a population × predator-cue design, so you can measure both evolved differences and plasticity.

## Main questions

1. **Does plasticity help or hinder adaptive evolution?** Theory predicts both. The authors tested it in a natural transplant experiment.
2. **Did brain gene expression evolve quickly?** They compared the introduced populations with their source after only 3–4 generations in cichlid-free streams.
3. **Does the initial plastic response point the same way as later evolution?** For example, when Guanapo fish are raised without predators, do their expression changes go in the same direction as the evolved changes, or the opposite?

## Key contrasts

- `gh` vs `ul`/`ll`: evolved change after introduction
- `S1` vs `S2` within `gh`: ancestral plasticity
- `tl`: reference for the adapted low-predation state
