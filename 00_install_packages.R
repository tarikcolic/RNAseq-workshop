###############################################################################
#
#   TRANSCRIPTOMICS SUMMER SCHOOL — SARAJEVO
#   RUN THIS BEFORE THE WORKSHOP
#
#   ---------------------------------------------------------------------------
#   DO THIS AT HOME, NOT ON THE FIRST MORNING. It takes 30-60 minutes on a
#   decent connection and considerably longer on conference wifi. If it
#   fails, you want to find out with a week to spare, not at 09:00 on Monday.
#   ---------------------------------------------------------------------------
#
#   BEFORE YOU START
#   ----------------
#   1. Install R, version 4.3 or newer:  https://cran.r-project.org
#      Check with:  R.version.string
#   2. Install RStudio Desktop:
#      https://posit.co/download/rstudio-desktop/
#   3. macOS users: also install the Xcode command line tools, or some
#      packages will fail to compile. In Terminal:
#           xcode-select --install
#      Windows users: install Rtools for your R version from
#           https://cran.r-project.org/bin/windows/Rtools/
#
#   HOW TO RUN
#   ----------
#   Open this file in RStudio and click "Source" (or Ctrl/Cmd + Shift + S).
#   Answer "no" if asked whether to install from source when a binary is
#   available. Answer "yes" if asked to create a personal library.
#
#   Then run scripts/00_check_installation.R, which prints a pass/fail
#   table. Send your instructor a screenshot of that table.
#
###############################################################################


# ---- 0. Basic checks --------------------------------------------------------

R.version.string

if (getRversion() < "4.3.0")
  stop("Your R is too old (", getRversion(), "). ",
       "Please install R >= 4.3 from https://cran.r-project.org and try again.")

# Where will packages be installed? If this is a system directory you may
# hit permission errors; R will normally offer to create a personal library
# instead, and you should accept.
.libPaths()

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Some of these downloads are large; R's default 60-second timeout is not
# enough on a slow connection.
options(timeout = max(1200, getOption("timeout")))

# Prefer pre-compiled binaries. Compiling Seurat's dependencies from source
# takes a very long time and is the most common cause of a failed install.
options(install.packages.compile.from.source = "never")


# ---- 1. A helper that installs only what is missing -------------------------
# Re-running the whole script should be cheap. This function skips anything
# already present, and reports clearly at the end.

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) {
    message("  all present: ", paste(pkgs, collapse = ", "))
    return(invisible(character(0)))
  }
  message("  installing: ", paste(missing, collapse = ", "))
  if (bioc) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing)
  }
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still))
    warning("FAILED to install: ", paste(still, collapse = ", "), call. = FALSE)
  invisible(still)
}

failures <- character(0)


# ---- 2. BiocManager, the gateway to Bioconductor ---------------------------
# CRAN hosts general-purpose R packages. Bioconductor is a separate,
# release-synchronised repository for biological data analysis, and it has
# its own installer. DESeq2, clusterProfiler and the annotation packages all
# live there.

message("\n== BiocManager ==")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)
BiocManager::version()

# Bioconductor versions are tied to R versions. If this warns about a
# mismatch, the R upgrade in step 0 is the fix.
BiocManager::valid(fix = FALSE, silent = TRUE)


# ---- 3. CRAN: general tools, tidyverse, plotting ---------------------------
message("\n== CRAN: core tools ==")
failures <- c(failures, install_if_missing(c(
  "tidyverse",      # dplyr, ggplot2, tidyr, readr, tibble, stringr, purrr
  "data.table",     # fast file reading, used for big count tables
  "R.utils",        # gzip handling for some downloads
  "remotes"         # installing from GitHub, if needed later
)))

message("\n== CRAN: figures ==")
failures <- c(failures, install_if_missing(c(
  "pheatmap",       # the simple, reliable heatmap (Day 3, Day 5)
  "RColorBrewer",   # colour-blind-safe palettes
  "viridis",        # perceptually uniform, colour-blind safe
  "ggrepel",        # non-overlapping text labels (volcano plots)
  "patchwork",      # compose multi-panel figures with + / | and /
  "cowplot",        # alternative composer; Seurat uses it internally
  "gridExtra"       # grid.arrange, used in one Day 3 figure
)))


# ---- 4. Bioconductor: bulk RNA-seq ----------------------------------------
message("\n== Bioconductor: differential expression (REQUIRED) ==")
failures <- c(failures, install_if_missing(c(
  "DESeq2",                 # the main engine for Days 2-3
  "apeglm",                 # log2 fold change shrinkage (lfcShrink)
  "SummarizedExperiment"    # the container the example datasets come in
), bioc = TRUE))

# Used only by the optional extension sections. A failure here does not
# stop you doing the workshop, so we do not add these to `failures`.
message("\n== Bioconductor: differential expression (optional) ==")
install_if_missing(c(
  "ashr",                   # alternative shrinkage, for arbitrary contrasts
  "edgeR",                  # the main alternative to DESeq2 (Day 2, ext.)
  "limma",                  # linear models; removeBatchEffect (Day 3, ext.)
  "EnhancedVolcano",        # an alternative volcano-plot package
  "ComplexHeatmap"          # a more powerful heatmap package
), bioc = TRUE)


# ---- 5. Bioconductor: functional interpretation ---------------------------
message("\n== Bioconductor: enrichment and annotation ==")
failures <- c(failures, install_if_missing(c(
  "clusterProfiler",  # GO/KEGG over-representation and GSEA
  "enrichplot",       # dotplot, emapplot, cnetplot for the above
  "fgsea",            # fast GSEA
  "GO.db",            # the Gene Ontology itself (term names and structure)
  "AnnotationDbi",    # the interface to all the org.*.eg.db packages
  "org.Hs.eg.db",     # human gene annotation      (groups 1, 3, 5)
  "org.Dm.eg.db",     # Drosophila gene annotation (group 2)
  "org.Mm.eg.db"      # mouse gene annotation      (single-cell groups 3, 5)
), bioc = TRUE))

# msigdbr gives access to the MSigDB gene set collections. Recent versions
# download the collections on first use, which needs internet and a few
# minutes. It is OPTIONAL — the Day 3 script falls back to GO gene sets
# built from the annotation packages above, which work entirely offline.
message("\n== msigdbr (optional) ==")
install_if_missing("msigdbr")


# ---- 6. Bulk example datasets (Days 2-3) ----------------------------------
# These are Bioconductor DATA packages: each ships a real published count
# matrix plus its sample table. Together about 60 MB.
#
# You only strictly need YOUR group's dataset, but installing all of them
# means you can look at what the other groups are doing, and it makes the
# round-up discussions much more useful.
message("\n== example datasets: bulk ==")
failures <- c(failures, install_if_missing(c(
  "airway",          # group 1 and group 5: human airway, dexamethasone
  "pasilla",         # group 2: Drosophila, pasilla knockdown
  "parathyroidSE",   # group 3: human parathyroid, DPN/OHT, 24/48 h
  "fission"          # group 4: S. pombe osmotic stress time course
), bioc = TRUE))


# ---- 7. Single cell (Days 4-5) --------------------------------------------
# Seurat has many dependencies. This is the slowest part of the install;
# expect 10-25 minutes. If it fails, the cause is almost always a missing
# compiler (see step 0) — fix that and re-run.
message("\n== Seurat (REQUIRED) ==")
failures <- c(failures, install_if_missing(c(
  "Seurat",          # the single-cell analysis framework
  "Matrix"           # sparse matrices; used directly, not just via Seurat
)))

message("\n== single cell, optional extensions ==")
install_if_missing("clustree")           # clustering trees (Day 5, section 2.2)
install_if_missing(c(
  "SingleCellExperiment",  # the Bioconductor single-cell container
  "scDblFinder",           # doublet detection (Day 4, section 9.3)
  "SingleR",               # automated annotation (Day 5, section 9.1)
  "celldex"                # reference datasets for SingleR
), bioc = TRUE)
install_if_missing("harmony")            # integration (Day 5, section 9.4)


# ---- 8. Report --------------------------------------------------------------
message("\n", strrep("=", 70))
if (length(failures) == 0) {
  message("All required packages installed successfully.")
} else {
  message("The following packages FAILED to install:\n  ",
          paste(unique(failures), collapse = ", "),
          "\n\nWhat to do:\n",
          "  1. Check step 0: do you have Rtools (Windows) or the Xcode\n",
          "     command line tools (macOS) installed?\n",
          "  2. Re-run this script — transient download failures are common\n",
          "     and re-running only retries what is missing.\n",
          "  3. Read the FIRST error in the console output, not the last.\n",
          "  4. If it still fails, email your instructor the output of\n",
          "     sessionInfo() and the first error message.")
}
message(strrep("=", 70))

message("\nNow run scripts/00_check_installation.R and send your instructor",
        "\na screenshot of the table it prints.\n")

sessionInfo()
