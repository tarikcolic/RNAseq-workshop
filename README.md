# RNAseq Workshop

This repository contains the materials for a hands-on transcriptomics workshop focused on RNA-seq and single-cell RNA-seq analysis in R/RStudio. The exercises are designed to take participants from foundational R skills through differential expression analysis, biological interpretation, and visualization of high-dimensional genomic data.

## Workshop leadership and credit

This workshop was led by Nidal Karagic, PhD., a Principal Investigator based in Helsinki. His guidance shaped the scientific direction of the course, the structure of the practical exercises, and the interpretation of the biological questions addressed throughout the week. We are grateful for his leadership and mentorship.

## Learning goals

By the end of the workshop, participants were expected to:

- work confidently in R and RStudio
- import, explore, and tidy count and metadata tables
- perform quality control and preprocessing for bulk RNA-seq data
- run differential gene expression analysis with DESeq2
- visualize transcriptomic patterns with PCA, heatmaps, and volcano plots
- interpret pathway and gene-set enrichment results
- analyze single-cell RNA-seq data using Seurat
- apply the workflow to real biological datasets and present findings

## Repository structure

- `00_install_packages.R` — package installation script for the required R/Bioconductor dependencies
- `00_check_installation.R` — environment verification script
- `01_R_basics.R` — foundational R exercises
- `Day1/` — Day 1 materials and practice scripts
- `Day2/` — bulk RNA-seq dataset and workshop lecture material
- `Day3/` — visualization and gene ontology enrichment analysis
- `Day4/` — single-cell RNA-seq Seurat workflow
- `Day5/` — final workshop lecture and feedback materials
- `Data/` — example datasets used throughout the course
- `analyses/` — output plots and analysis results
- `output/` — processed materials and exported tables
- `Group_project/` — group exercise scripts and project work

## What was learned each day

### Day 1: Foundations in R and reproducible analysis

The first day introduced the R environment and the core habits needed for computational biology. Participants learned how to:

- set up and navigate an RStudio project
- use variables, vectors, and data structures
- import tabular data and inspect sample metadata
- perform calculations and summarise numeric information
- generate plots to explore patterns in data
- write reproducible scripts and work interactively in the console

This day focused on building confidence with the scripting language that underpins the rest of the workshop.

### Day 2: Bulk RNA-seq analysis with DESeq2

The second day introduced the bulk RNA-seq workflow from count tables to differential expression testing. Participants learned how to:

- load expression count matrices and sample annotation tables
- understand experimental design and model formulae
- normalize counts and filter low-information genes
- fit DESeq2 models for biological comparisons
- interpret differential expression results including log2 fold change and adjusted p-values
- identify genes that are significantly up- or down-regulated across conditions

This is the core statistical foundation for DE analysis in transcriptomics.

### Day 3: Visualization and functional interpretation

The third day focused on showing how to make results understandable and biologically interpretable. Participants learned how to:

- perform variance-stabilizing transformation and assess sample structure
- create PCA plots and sample-distance heatmaps
- make volcano plots and expression heatmaps
- translate gene identifiers between annotation systems
- run GO enrichment and GSEA analyses
- interpret enriched biological processes and pathways

The emphasis here was not only on generating plots, but also on choosing the right visualization to support biological conclusions.

### Day 4: Single-cell RNA-seq with Seurat

The fourth day moved into single-cell transcriptomics. Participants learned the standard Seurat workflow, including:

- loading 10x Genomics count matrices
- performing quality control for cells and genes
- normalizing and scaling data
- identifying variable features
- running PCA and selecting informative dimensions
- clustering cells and visualizing them with UMAP
- identifying marker genes for each cluster
- interpreting cell type structure in a scRNA-seq dataset

This day gave participants a practical introduction to cellular heterogeneity and cluster-based analysis.

### Day 5: Integration, interpretation, and project work

The final day brought the workflow together through project-based analysis and synthesis. Participants applied the skills learned earlier to their own datasets and discussed:

- how to choose appropriate analysis strategies for a given biological question
- how to interpret pre-processing, clustering, and differential expression results critically
- how to connect statistical output to biological meaning
- how to present transcriptomics findings clearly in a research context

This day focused on integration, interpretation, and communication of results.

## Suggested workflow

To work through the examples in this repository:

1. Open the project in RStudio using `RNAseq_workshop.Rproj`
2. Run `00_install_packages.R` before the workshop begins
3. Check the environment with `00_check_installation.R`
4. Follow the scripts in order, starting with `01_R_basics.R` and then the Day-specific exercises
5. Use the `Data/` and `analyses/` folders to connect the code to the datasets and results

## Data sources

The repository includes multiple example datasets used for teaching, including bulk RNA-seq data and a single-cell PBMC dataset. These examples were chosen to illustrate common analysis tasks and real-world biological interpretation scenarios.

## Acknowledgements

This workshop brought together practical bioinformatics training, R-based analysis, and project-based learning in transcriptomics. Special thanks to Nidal Karagic, PhD., whose leadership and scientific mentorship guided the workshop and inspired the participants throughout the week.

## License

This repository is intended for educational use within the workshop context. Please credit the original authors and instructors when reusing the material.

