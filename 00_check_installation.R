###############################################################################
#
#   TRANSCRIPTOMICS SUMMER SCHOOL — SARAJEVO
#   Installation check
#
#   Run this AFTER scripts/00_install_packages.R.
#   Open it in RStudio and click "Source".
#
#   It prints a table with one row per requirement. Send your instructor a
#   screenshot of the whole console output. If anything says FAIL, the
#   message at the bottom tells you what to do.
#
###############################################################################

cat("\n")
cat(strrep("=", 74), "\n")
cat("  TRANSCRIPTOMICS SUMMER SCHOOL — INSTALLATION CHECK\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat(strrep("=", 74), "\n\n")


# ---- 1. R version -----------------------------------------------------------
r_ok <- getRversion() >= "4.3.0"
cat(sprintf("%-14s %-30s %s\n", "R version", as.character(getRversion()),
            if (r_ok) "OK" else "FAIL — need >= 4.3.0"))
cat(sprintf("%-14s %-30s\n", "platform", R.version$platform))
cat("\n")


# ---- 2. Packages ------------------------------------------------------------
# Grouped by the day that needs them, so a failure tells you what you will
# be unable to do.

# REQUIRED = the scripts will not run without it.
groups <- list(
  "Day 1 — R basics" = c("tidyverse", "dplyr", "ggplot2", "readr", "tibble"),
  "Days 2-3 — bulk"  = c("DESeq2", "apeglm", "SummarizedExperiment"),
  "Days 2-3 — plots" = c("pheatmap", "RColorBrewer", "ggrepel", "patchwork",
                         "gridExtra"),
  "Day 3 — pathways" = c("clusterProfiler", "enrichplot", "fgsea", "GO.db",
                         "AnnotationDbi", "org.Hs.eg.db", "org.Dm.eg.db"),
  "Bulk datasets"    = c("airway", "pasilla", "parathyroidSE", "fission"),
  "Days 4-5 — single cell" = c("Seurat", "Matrix")
)

# OPTIONAL = used only by the optional extension sections at the end of each
# script, or by a fallback path. Missing ones are reported but do not fail
# the check. Install them if you can; the workshop works without them.
optional <- c("ashr", "limma", "edgeR",          # Day 2-3 extensions
              "clustree",                         # Day 5 clustering tree
              "SingleCellExperiment", "scDblFinder", "SingleR", "celldex",
              "harmony",                          # Day 4-5 extensions
              "EnhancedVolcano", "ComplexHeatmap", "org.Mm.eg.db",
              "msigdbr", "cowplot", "viridis", "data.table")

n_fail <- 0L

for (g in names(groups)) {
  cat(g, "\n", strrep("-", nchar(g)), "\n", sep = "")
  for (p in groups[[g]]) {
    ok <- requireNamespace(p, quietly = TRUE)
    if (!ok) n_fail <- n_fail + 1L
    v <- if (ok) as.character(utils::packageVersion(p)) else "-"
    cat(sprintf("  %-24s %-12s %s\n", p, v, if (ok) "OK" else "FAIL"))
  }
  cat("\n")
}

cat("Optional (nice to have, not required)\n")
cat(strrep("-", 37), "\n", sep = "")
for (p in optional) {
  ok <- requireNamespace(p, quietly = TRUE)
  v  <- if (ok) as.character(utils::packageVersion(p)) else "-"
  cat(sprintf("  %-24s %-12s %s\n", p, v, if (ok) "ok" else "missing"))
}
cat("\n")


# ---- 3. Do the workshop's own helper files load? ---------------------------
# This also confirms that you have the folder in one piece and that R can
# find the workshop root from wherever your session is sitting.
cat("Workshop files\n", strrep("-", 14), "\n", sep = "")

root_found <- tryCatch({
  p <- normalizePath(getwd(), winslash = "/")
  hit <- NA_character_
  for (i in 0:5) {
    if (file.exists(file.path(p, "scripts", "workshop_data.R"))) { hit <- p; break }
    if (identical(dirname(p), p)) break
    p <- dirname(p)
  }
  hit
}, error = function(e) NA_character_)

if (is.na(root_found)) {
  cat("  workshop root            -            FAIL\n")
  cat("    -> Open 24_Sarajevo_Summer_School.Rproj, or setwd() into the\n")
  cat("       workshop folder, then run this script again.\n")
  n_fail <- n_fail + 1L
} else {
  cat(sprintf("  %-24s %s\n", "workshop root", root_found))
  for (f in c("scripts/workshop_data.R", "scripts/workshop_sc_data.R")) {
    ex <- file.exists(file.path(root_found, f))
    cat(sprintf("  %-24s %-12s %s\n", basename(f), "", if (ex) "OK" else "FAIL"))
    if (!ex) n_fail <- n_fail + 1L
  }
}
cat("\n")


# ---- 4. A real smoke test --------------------------------------------------
# Installed is not the same as working. Load a real dataset and fit a real
# (tiny) model. This catches broken compiled dependencies, which "is the
# package installed" checks cannot.
cat("Smoke test — can we actually run an analysis?\n")
cat(strrep("-", 44), "\n", sep = "")

smoke <- tryCatch({
  suppressPackageStartupMessages({
    library(DESeq2); library(airway)
  })
  data("airway", package = "airway")
  se <- airway
  cts <- SummarizedExperiment::assay(se)
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  cd$dex <- relevel(factor(cd$dex), ref = "untrt")
  dds <- DESeqDataSetFromMatrix(cts, cd, ~ dex)
  # keep a small, well-expressed subset so this takes seconds, not minutes
  dds <- dds[rowSums(counts(dds) >= 50) == ncol(dds), ][1:400, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("dex", "trt", "untrt"))
  sprintf("DESeq2 fitted %d genes on %d samples; %d at padj<0.1",
          nrow(res), ncol(dds), sum(res$padj < 0.1, na.rm = TRUE))
}, error = function(e) paste("FAILED:", conditionMessage(e)))

cat("  bulk (DESeq2 + airway)   ", smoke, "\n")
if (grepl("^FAILED", smoke)) n_fail <- n_fail + 1L

smoke_sc <- tryCatch({
  suppressPackageStartupMessages(library(Seurat))
  set.seed(1)
  m <- matrix(rpois(200 * 120, 1.2), nrow = 200,
              dimnames = list(paste0("G", 1:200), paste0("C", 1:120)))
  o <- CreateSeuratObject(m)
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, nfeatures = 100, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = 10, verbose = FALSE)
  o <- FindNeighbors(o, dims = 1:5, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.5, verbose = FALSE)
  sprintf("Seurat v%s ran the full pipeline on %d synthetic cells (%d clusters)",
          as.character(packageVersion("Seurat")), ncol(o),
          length(unique(Idents(o))))
}, error = function(e) paste("FAILED:", conditionMessage(e)))

cat("  single cell (Seurat)     ", smoke_sc, "\n")
if (grepl("^FAILED", smoke_sc)) n_fail <- n_fail + 1L

# Can we write files? Needed for every figure and results table.
write_ok <- tryCatch({
  tf <- file.path(tempdir(), "workshop_write_test.csv")
  write.csv(data.frame(a = 1:2), tf, row.names = FALSE)
  ok <- file.exists(tf); unlink(tf); ok
}, error = function(e) FALSE)
cat("  can write files          ", if (write_ok) "OK" else "FAIL", "\n")
if (!write_ok) n_fail <- n_fail + 1L

# Can we reach the internet? Needed once, on Day 4, to fetch the
# single-cell data. Not fatal now — the instructor has the files.
net_ok <- tryCatch({
  con <- url("https://cloud.r-project.org", open = "rb")
  on.exit(close(con), add = TRUE)
  length(readBin(con, "raw", 64)) > 0
}, error = function(e) FALSE)
cat("  internet reachable       ", if (net_ok) "OK" else
    "no (fine for now; needed on Day 4)", "\n\n")


# ---- 5. Verdict -------------------------------------------------------------
cat(strrep("=", 74), "\n")
if (n_fail == 0 && r_ok) {
  cat("  PASS — you are ready for the workshop. See you in Sarajevo.\n")
} else {
  cat(sprintf("  %d PROBLEM(S) FOUND.\n\n", n_fail + !r_ok))
  cat("  1. If R itself is too old, install R >= 4.3 first; everything else\n")
  cat("     follows from that.\n")
  cat("  2. Otherwise, re-run scripts/00_install_packages.R. It only retries\n")
  cat("     what is missing, and transient download failures are common.\n")
  cat("  3. If a package still fails, you are probably missing a compiler:\n")
  cat("       Windows: install Rtools\n")
  cat("       macOS:   run  xcode-select --install  in Terminal\n")
  cat("  4. Then email your instructor this console output together with\n")
  cat("     the sessionInfo() below. Do this BEFORE the workshop starts.\n")
}
cat(strrep("=", 74), "\n\n")

sessionInfo()
