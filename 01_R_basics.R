###############################################################################
#
#   TRANSCRIPTOMICS SUMMER SCHOOL — SARAJEVO
#   Day 1, afternoon:  R as a research tool
#
#   GOAL OF THIS SESSION
#   --------------------
#   By the end of today you should be able to load a table, look at it,
#   filter it, compute summaries from it, and plot it — without needing to
#   look up the syntax. That is all the R you need for the rest of the week.
#
#   HOW TO USE THIS FILE
#   --------------------
#   Do NOT run the whole script at once. Work through it line by line.
#   In RStudio: put the cursor on a line and press
#        Ctrl + Enter        (Windows/Linux)
#        Cmd  + Enter        (macOS)
#   The line is sent to the console and executed. Read the output BEFORE
#   moving on. If you do not understand the output, that is the moment to
#   ask — not three sections later.
#
#   Lines starting with '#' are comments: R ignores them. They are for you.
#   Sections are marked with  # ---- NAME ----  which RStudio turns into a
#   navigable outline (bottom-left of the editor pane, or Ctrl/Cmd + Shift + O).
#
#   Every section ends with an EXERCISE. Do them. Reading code is not the
#   same as writing it, and only one of the two transfers to tomorrow.
#
###############################################################################


# ---- 0. Set-up and orientation ----------------------------------------------

# Where am I? R always has a "working directory": the folder it looks in
# when you give it a file name without a full path.
getwd()

# If you opened this file via the .Rproj file (recommended), the working
# directory is already the workshop folder and you should never need setwd().
# Avoid setwd() in scripts you send to other people: the path only exists
# on your machine.

# What is loaded right now? (empty at the start of a session)
ls()

# R's four panes in RStudio:
#   Source (this file) | Console (where code runs)
#   Environment (your objects) | Files / Plots / Help / Packages

# The single most useful command in R:
?mean          # opens the help page for mean(). Scroll to "Examples".

# Keep this habit: when a function surprises you, read its help page before
# guessing. Help pages tell you the arguments, the defaults, and what the
# function returns.



# ---- 1. R as a calculator; assignment ---------------------------------------

# R evaluates expressions and prints the result.
2 + 2
10 / 4
2^10
sqrt(144)
log(100)          # NOTE: log() is the NATURAL log (base e), not base 10
log10(100)        # base 10
log2(8)           # base 2 — you will use this constantly (fold changes)
exp(1)

# Integer division and remainder — occasionally handy
17 %/% 5          # how many whole 5s fit into 17
17 %%  5          # what is left over


## ASSIGNMENT -----------------------------------------------------------------
# To keep a value, give it a name. The assignment operator is  <-
# (You may also see '='. Use '<-' for assignment; it is the R convention
#  and avoids ambiguity with function arguments.)

x <- 5
y <- 3

x                 # typing the name prints the object
x + y
z <- x * y        # assignment does NOT print; the value is stored silently
z

# Naming rules and conventions:
#   - start with a letter
#   - letters, numbers, '.' and '_' allowed
#   - R is CASE SENSITIVE:  gene, Gene and GENE are three different objects
#   - do not name objects after existing functions (c, mean, data, df, T, F)
#     It is legal but it will eventually confuse you badly.

# Good names describe content:
mean_expression <- 42.7
n_samples       <- 8

# Overwrite freely — R does not warn you. This is a common source of bugs.
x <- 100
x


## THE CLASSES OF SINGLE VALUES ----------------------------------------------
# Every object has a class. Ask with class().

class(5)                  # "numeric"   — any number
class(5L)                 # "integer"   — the L forces an integer
class("hello")            # "character" — text, ALWAYS in quotes
class(TRUE)               # "logical"   — TRUE / FALSE (capitals, no quotes)
class(factor("control"))  # "factor"    — categorical; more in section 5

# Special values you must be able to recognise:
NA        # Not Available — a MISSING value. Extremely common in real data.
NULL      # the empty object — "nothing is here"
NaN       # Not a Number   — e.g. 0/0
Inf       # infinity       — e.g. 1/0

0/0
1/0
is.na(NA)                 # test for missingness — never use  == NA


## EXERCISE 1 ----------------------------------------------------------------
# a) Store the number of samples in your future experiment in `n`.
# b) Store the sequencing depth per sample (say 25 million) in `depth`.
# c) Compute the total number of reads and store it in `total_reads`.
# d) Express `total_reads` in billions, rounded to two decimals.
#    Hint: ?round



# ---- 2. Vectors and vectorised operations -----------------------------------

# A VECTOR is an ordered collection of values OF THE SAME TYPE.
# It is the fundamental data structure in R. c() = "combine".

counts <- c(120, 340, 55, 890, 12)
counts

genes  <- c("TP53", "MYC", "GAPDH", "ACTB", "EGFR")
genes

length(counts)            # how many elements
class(counts)
class(genes)

# Handy ways to build vectors:
1:10                      # integer sequence
seq(0, 100, by = 25)      # explicit steps
seq(0, 1, length.out = 5) # explicit number of elements
rep("ctrl", 3)            # repeat a value
rep(c("ctrl", "trt"), each = 3)   # <- exactly how you build a condition vector
rep(c("ctrl", "trt"), times = 3)  # note the difference!


## VECTORISATION -------------------------------------------------------------
# Operations apply to EVERY element at once. There is no loop, and you
# should not write one. This is the single most important idea in R.

counts * 2
counts + 1000
log2(counts + 1)          # the +1 is a "pseudocount": log2(0) is -Inf

# Two vectors of the same length operate element by element:
ctrl <- c(100, 200, 300, 400)
trt  <- c(150, 180, 600, 400)

trt - ctrl                # difference per gene
trt / ctrl                # ratio per gene
log2(trt / ctrl)          # log2 fold change — THE standard expression metric

# Why log2 fold change?
#   - it is symmetric: 2x up = +1, 2x down = -1
#     (whereas on the ratio scale, up is 2 and down is 0.5 — not symmetric)
#   - it makes the distribution roughly symmetric around zero
#   - biologists can read it: +3 means 8-fold up


## CAUTION: RECYCLING -------------------------------------------------------
# If lengths differ, R silently REPEATS the shorter vector. This is
# occasionally useful and frequently a bug.
c(1, 2, 3, 4) + c(10, 20)       # becomes 1+10, 2+20, 3+10, 4+20
# c(1,2,3) + c(10,20)           # uncomment: warns, because 3 is not a multiple of 2

# Defensive habit: check lengths match before combining vectors.
length(ctrl) == length(trt)


## SUMMARY FUNCTIONS --------------------------------------------------------
sum(counts)
mean(counts)
median(counts)
sd(counts)                # standard deviation
var(counts)
min(counts); max(counts)
range(counts)
summary(counts)           # five-number summary + mean

# Missing values BLOCK most summaries by default — this is a feature,
# because it forces you to notice them:
counts_na <- c(120, 340, NA, 890, 12)
mean(counts_na)                    # NA
mean(counts_na, na.rm = TRUE)      # remove NAs and compute
sum(is.na(counts_na))              # how many are missing?


## NAMED VECTORS ------------------------------------------------------------
# Elements can carry names. Very useful: a count vector labelled by gene.
names(counts) <- genes
counts
counts["MYC"]             # access by name rather than position
sort(counts, decreasing = TRUE)


## EXERCISE 2 ---------------------------------------------------------------
# a) Create a vector `expr` with the values 5, 100, 2500, 0, 48, 7.
# b) Name its elements GENE1 ... GENE6.
# c) Compute the log2 of expr + 1.
# d) Which gene has the highest value? (Hint: ?which.max, and use the names)
# e) A second sample gives 8, 90, 5000, 3, 40, 14. Compute the log2 fold
#    change of sample 2 relative to sample 1 for every gene.



# ---- 3. Data frames; reading and writing files ------------------------------

# A DATA FRAME is a table: columns may have different types, but every
# column must have the same length. This is how metadata / sample sheets
# and results tables are represented.

samples <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4", "S5", "S6"),
  condition = c("control", "control", "control", "treated", "treated", "treated"),
  batch     = c("A", "B", "A", "B", "A", "B"),
  rin       = c(9.1, 8.7, 9.4, 7.9, 8.8, 9.0),
  libsize   = c(24.1e6, 19.8e6, 27.3e6, 22.0e6, 25.5e6, 18.2e6),
  stringsAsFactors = FALSE      # keep text as text; we will convert
)                               # to factors deliberately in section 5

samples


## INSPECT BEFORE YOU DO ANYTHING ELSE --------------------------------------
# These five calls should be automatic every time you load data.

dim(samples)        # rows, columns
nrow(samples); ncol(samples)
names(samples)      # column names (colnames() also works)
str(samples)        # STRUCTURE: types + first values. Use this constantly.
head(samples, 3)    # first rows
tail(samples, 2)    # last rows
summary(samples)    # per-column summary


## ACCESSING COLUMNS --------------------------------------------------------
samples$condition           # $ returns the column as a vector
samples[["condition"]]      # equivalent
samples[, "condition"]      # equivalent (base-R matrix-style)

# A column is just a vector, so everything from section 2 applies:
mean(samples$rin)
table(samples$condition)    # counts per category — extremely useful
table(samples$condition, samples$batch)   # cross-tabulation
                            # ^ THIS is how you check your design is balanced


## ADDING AND MODIFYING COLUMNS --------------------------------------------
samples$libsize_millions <- samples$libsize / 1e6
samples$low_rin          <- samples$rin < 8.5     # logical column
samples


## WRITING AND READING FILES -----------------------------------------------
# We write into a subfolder so we do not clutter the project root.
#
# NOTE where these files land: "output/..." is relative to the WORKING
# DIRECTORY, which is whatever getwd() said in section 0 — probably the
# workshop root if you opened the .Rproj. That is the whole point of
# relative paths, and also the most common source of "where did my file
# go?". If you cannot find the file, run getwd() and look there.
# (From Day 2 onward the scripts use a helper that resolves absolute paths,
#  so your results always land in the right day's folder.)
dir.create("output", showWarnings = FALSE)  # showWarnings=FALSE: fine if it exists

# Comma-separated values. row.names = FALSE avoids an unnamed first column,
# which is the most common source of "why is my table shifted by one".
write.csv(samples, file = "output/sample_metadata.csv", row.names = FALSE)

# Tab-separated is often safer in genomics (gene descriptions contain commas)
write.table(samples, file = "output/sample_metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Read it back. Arguments worth knowing:
#   header    — is the first line column names?
#   sep       — the delimiter
#   row.names — which column holds row names (often 1 for count matrices)
back <- read.csv("output/sample_metadata.csv", header = TRUE)
str(back)

# Reading a count matrix typically looks like this (gene IDs in column 1):
#   counts <- read.delim("counts.tsv", row.names = 1, check.names = FALSE)
#   counts <- as.matrix(counts)
# check.names = FALSE stops R from mangling sample names like "S1-treated"
# into "S1.treated".


## EXERCISE 3 ---------------------------------------------------------------
# a) Add a column `group_batch` that pastes condition and batch together,
#    separated by "_".  Hint: ?paste
# b) How many samples are in each condition x batch combination?
# c) Is this design balanced with respect to batch? Say why in a comment.
# d) Write the modified table to output/exercise3.tsv as tab-separated.



# ---- 4. Subsetting and logical indexing -------------------------------------

# This is the skill you will use most. Take it slowly.

## VECTORS: by position ----------------------------------------------------
counts                  # (still named, from section 2)
counts[1]               # first element. NOTE: R counts from 1, not 0.
counts[c(1, 3)]         # first and third
counts[2:4]             # a range
counts[-1]              # everything EXCEPT the first (negative = drop)


## VECTORS: by logical condition -------------------------------------------
# A comparison produces a logical vector of TRUE/FALSE ...
counts > 100
# ... and a logical vector can be used as an index: keep where TRUE.
counts[counts > 100]

# The comparison operators:
#   ==  equal          !=  not equal
#   <   less           <=  less or equal
#   >   greater        >=  greater or equal
#   %in% is a member of

genes %in% c("MYC", "EGFR")
counts[genes %in% c("MYC", "EGFR")]

# Combine conditions:
#   &   AND      |   OR      !   NOT
counts[counts > 50 & counts < 500]
counts[counts < 50 | counts > 500]

# NEVER use  == NA . Use is.na().
#   NA == NA  is NA, not TRUE — because "unknown equals unknown" is unknown.
NA == NA
is.na(NA)


## DATA FRAMES: [rows, columns] --------------------------------------------
# The comma is essential. Before it: rows. After it: columns.
# Leaving one side blank means "all of them".

samples[1, ]                      # first ROW, all columns
samples[, 2]                      # all rows, second COLUMN
samples[1:3, c("sample_id", "rin")]   # rows 1-3, two named columns

# Filtering rows by a condition on a column — the workhorse pattern:
samples[samples$condition == "treated", ]
samples[samples$rin > 8.8, ]
samples[samples$condition == "control" & samples$batch == "A", ]

# Just the sample IDs that pass:
samples$sample_id[samples$rin > 8.8]

# Sorting a data frame: order() returns the ORDER of the indices,
# which you then use to subset the rows.
samples[order(samples$rin), ]
samples[order(-samples$rin), ]                    # descending
samples[order(samples$condition, samples$rin), ]  # sort by two keys


## THE ORDERING TRAP -------------------------------------------------------
# This is the mistake that has caused real papers to be retracted.
# Your count matrix columns and your metadata rows MUST correspond.

count_cols <- c("S3", "S1", "S2", "S6", "S4", "S5")   # as they came from the pipeline
count_cols
samples$sample_id                                      # as in your sample sheet

# Are they in the same order? No.
all(count_cols == samples$sample_id)

# Reorder the metadata to MATCH the matrix (never the other way round,
# because the matrix columns are the ground truth of what was sequenced):
samples_ordered <- samples[match(count_cols, samples$sample_id), ]
samples_ordered
all(count_cols == samples_ordered$sample_id)   # must be TRUE before you continue

# Make this a hard check in every analysis script you ever write:
stopifnot(identical(count_cols, samples_ordered$sample_id))
# stopifnot() throws an error and stops the script if the condition is FALSE.
# It is the cheapest insurance in data analysis.


## EXERCISE 4 ---------------------------------------------------------------
# a) Extract all rows where libsize is below 23 million.
# b) Extract the rin values of the treated samples only.
# c) How many control samples have rin above 9.0?
# d) Sort the table by libsize, largest first, and show only sample_id
#    and libsize_millions.



# ---- 5. Factors — categorical variables (and why they bite) -----------------

# A factor stores a categorical variable as integer codes plus a set of
# LEVELS (the allowed categories). Every statistical model in R —
# including DESeq2 tomorrow — treats factors specially.

cond <- factor(c("control", "control", "treated", "treated"))
cond
levels(cond)         # "control" "treated"  — alphabetical by default
as.integer(cond)     # the underlying codes

# WHY THIS MATTERS:
# The FIRST level is the REFERENCE level. In a model, every other level is
# compared AGAINST the reference. So the first level decides the direction
# and meaning of your fold changes.

# Alphabetical order is rarely what you want biologically:
cond2 <- factor(c("untreated", "untreated", "drug", "drug"))
levels(cond2)        # "drug" first, because d < u alphabetically!
# With this coding, DESeq2 would report "untreated vs drug" — the fold
# changes would be inverted relative to how you would write the paper.

# Fix it in one of two ways:
cond2 <- relevel(cond2, ref = "untreated")
levels(cond2)

# or set the levels explicitly when you create the factor (preferred,
# because it is self-documenting):
cond3 <- factor(c("untreated", "untreated", "drug", "drug"),
                levels = c("untreated", "drug"))
levels(cond3)

# Dropping unused levels after filtering — a classic source of confusion:
sub <- cond3[cond3 == "drug"]
levels(sub)          # still has "untreated" as a level, with 0 observations!
levels(droplevels(sub))

# Converting a factor of numbers back to numbers: NEVER do as.numeric(f)
# directly — you get the integer CODES, not the values.
f <- factor(c("10", "20", "30"))
as.numeric(f)                     # 1 2 3   <-- WRONG
as.numeric(as.character(f))       # 10 20 30 <-- correct


## EXERCISE 5 ---------------------------------------------------------------
# a) Convert samples$condition into a factor with "control" as the
#    reference level, and store it back in the data frame.
# b) Confirm with levels().
# c) Convert samples$batch into a factor. How many levels does it have?



# ---- 6. Writing your own functions ------------------------------------------

# Whenever you copy-paste code and change one thing, write a function
# instead. Functions are how you stop making transcription errors.

# ANATOMY:
#   name <- function(arguments) { body; return value }
# The value of the last expression is returned automatically, but writing
# return() explicitly is clearer for beginners.

cpm <- function(counts_vector) {
  # Counts Per Million: the simplest normalisation for library size.
  # counts_vector : numeric vector of raw counts for ONE sample
  # returns       : numeric vector of CPM values
  total <- sum(counts_vector)
  return(counts_vector / total * 1e6)
}

raw <- c(5, 200, 1000, 0, 45)
cpm(raw)
sum(cpm(raw))            # sums to 1e6 by construction — a good sanity check


# DEFAULT ARGUMENTS: give sensible defaults so the common case is short.
log_transform <- function(x, pseudocount = 1, base = 2) {
  # log-transform expression values after adding a pseudocount, so that
  # zeros do not become -Inf.
  log(x + pseudocount, base = base)
}

log_transform(raw)                    # uses the defaults
log_transform(raw, pseudocount = 0.5)
log_transform(raw, base = 10)


# DEFENSIVE FUNCTIONS: check your inputs. Errors close to the mistake are
# far cheaper to debug than errors three steps later.
log2_fc <- function(treated, control) {
  if (length(treated) != length(control)) {
    stop("treated and control must have the same length")
  }
  if (any(treated < 0, na.rm = TRUE) || any(control < 0, na.rm = TRUE)) {
    warning("negative values found — are these really counts?")
  }
  log2((treated + 1) / (control + 1))
}

log2_fc(c(100, 400, 50), c(100, 100, 200))
# log2_fc(c(1,2), c(1,2,3))     # uncomment to see the error


## SCOPE — a short but important warning ----------------------------------
# Objects created inside a function do NOT exist outside it.
demo <- function() {
  hidden <- 99
  hidden
}
demo()
# hidden        # uncomment: Error — object 'hidden' not found. This is good:
#               # it means functions cannot accidentally clobber your data.


## EXERCISE 6 ---------------------------------------------------------------
# a) Write a function `n_detected(x, min_count = 1)` that returns how many
#    elements of x are greater than or equal to min_count.
# b) Test it on c(0, 0, 3, 17, 0, 250) with the default and with min_count = 10.
# c) Write a function `summarise_gene(x)` that returns a NAMED vector with
#    the mean, median and sd of x. Hint: c(mean = ..., median = ..., sd = ...)



# ---- 7. Iteration: apply, sapply, and loops ---------------------------------

# Genomics data are matrices, and you very often want a per-gene (row) or
# per-sample (column) summary. Do not write a for-loop for this.

# Build a small fake count matrix: 6 genes x 4 samples
set.seed(42)          # makes the random numbers reproducible — ALWAYS do this
m <- matrix(rpois(24, lambda = 50), nrow = 6,
            dimnames = list(paste0("GENE", 1:6), paste0("S", 1:4)))
m

# Fast built-ins for the common cases:
rowSums(m)      # total counts per gene
colSums(m)      # library size per sample  <- you will use this tomorrow
rowMeans(m)
colMeans(m)

# apply(X, MARGIN, FUN):  MARGIN = 1 means rows, 2 means columns.
apply(m, 1, max)                      # maximum per gene
apply(m, 2, median)                   # median per sample
apply(m, 1, function(g) sum(g > 45))  # per gene: in how many samples > 45

# sapply() over a vector or list, returning a simplified result:
sapply(1:6, function(i) i^2)
sapply(as.data.frame(m), sd)          # sd per column

# lapply() returns a LIST (use when results are not all the same shape):
res <- lapply(as.data.frame(m), summary)
res$S1


## FOR LOOPS — when you do need them --------------------------------------
# Use a loop when each step has a SIDE EFFECT (writing a file, printing,
# saving a plot) rather than returning a value.
for (s in colnames(m)) {
  cat(s, ": library size =", sum(m[, s]), "\n")
}
# cat() prints without quotes; "\n" is a newline.

# If you loop to BUILD a result, pre-allocate — do not grow a vector.
out <- numeric(ncol(m))                # pre-allocated, correct length
names(out) <- colnames(m)
for (i in seq_along(colnames(m))) {    # seq_along is safer than 1:n
  out[i] <- sum(m[, i] > 45)
}
out
# ... but note that apply(m, 2, function(x) sum(x > 45)) does the same in
# one readable line. Prefer that.


## EXERCISE 7 ---------------------------------------------------------------
# a) Using apply(), compute the coefficient of variation (sd/mean) for each
#    gene in m.
# b) Which gene is most variable?
# c) Convert m to CPM using your cpm() function from section 6 applied
#    column-wise. Check that each column of the result sums to 1e6.



# ---- 8. dplyr: a readable grammar for tables --------------------------------

# Base R can do everything, but dplyr code reads like a description of what
# you want. In genomics you will mostly use it on results tables.

# install.packages("tidyverse")   # only once, if not already installed
library(dplyr)

# THE PIPE  |>  takes the thing on the left and passes it as the first
# argument to the function on the right.
#   x |> f()        is the same as   f(x)
#   x |> f(y)       is the same as   f(x, y)
# (The tidyverse pipe %>% behaves the same way; |> is now base R.)

# A more realistic object: a small differential expression results table.
set.seed(1)
res_tbl <- data.frame(
  gene      = paste0("GENE", sprintf("%03d", 1:200)),
  baseMean  = round(exp(rnorm(200, mean = 5, sd = 2)), 1),
  log2FC    = round(rnorm(200, mean = 0, sd = 1.3), 3),
  pvalue    = runif(200)^3                 # skewed towards small values
)
res_tbl$padj <- p.adjust(res_tbl$pvalue, method = "BH")  # FDR correction
head(res_tbl)

# THE SIX VERBS YOU NEED -------------------------------------------------

# filter()  — keep rows matching a condition
res_tbl |> filter(padj < 0.1) |> head()

# select()  — keep or drop columns
res_tbl |> select(gene, log2FC, padj) |> head()

# mutate()  — create or modify columns
res_tbl |>
  mutate(direction = ifelse(log2FC > 0, "up", "down")) |>
  head()

# arrange() — sort rows
res_tbl |> arrange(padj) |> head()
res_tbl |> arrange(desc(abs(log2FC))) |> head()

# group_by() + summarise() — split, compute, recombine
res_tbl |>
  mutate(direction = ifelse(log2FC > 0, "up", "down")) |>
  group_by(direction) |>
  summarise(n          = n(),
            mean_log2FC = mean(log2FC),
            n_signif    = sum(padj < 0.1))

# count() — the shortcut for group_by + tally
res_tbl |>
  mutate(signif = padj < 0.1) |>
  count(signif)


# CHAINING — this is the real payoff. Read it top to bottom as a sentence:
# "take the results, keep significant genes with a decent effect size,
#  add a direction label, sort by effect size, and show the top 10."
top_hits <- res_tbl |>
  filter(padj < 0.1, abs(log2FC) > 1) |>
  mutate(direction = ifelse(log2FC > 0, "up", "down"),
         fold      = round(2^abs(log2FC), 2)) |>
  arrange(desc(abs(log2FC))) |>
  select(gene, log2FC, fold, direction, padj) |>
  head(10)

top_hits

# The same thing in base R, for comparison — correct, but harder to read
# and much harder to modify:
tmp <- res_tbl[res_tbl$padj < 0.1 & abs(res_tbl$log2FC) > 1, ]
tmp <- tmp[order(-abs(tmp$log2FC)), ]
head(tmp[, c("gene", "log2FC", "padj")], 10)

# Note: dplyr operations RETURN a new table; they do not modify the original.
# If you want to keep the result, assign it (as we did with top_hits).


## EXERCISE 8 ---------------------------------------------------------------
# a) How many genes have padj < 0.05?
# b) Build a table of the 15 most significant genes showing gene, log2FC
#    and padj, sorted by padj.
# c) Add a column `expression_class` that is "low" if baseMean < 50,
#    "medium" if 50-500, and "high" if > 500. Hint: ?case_when
# d) Using group_by()/summarise(), report how many significant genes
#    (padj < 0.1) fall in each expression_class. What does the pattern
#    tell you about statistical power vs expression level?



# ---- 9. ggplot2: building a figure layer by layer ---------------------------

library(ggplot2)

# THE GRAMMAR
#   ggplot(data, aes(...))   the data, and which column maps to which
#                            visual property (aesthetic)
#   + geom_*()               what to draw
#   + scale_*()              how values become colours / positions / sizes
#   + labs()                 titles and axis labels
#   + theme_*()              overall appearance
# Layers are combined with '+'.


## 9.1 A minimal plot ----------------------------------------------------
ggplot(res_tbl, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point()

# Build it up. Note where each element appears in the output.
ggplot(res_tbl, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(alpha = 0.5, size = 1.8) +          # alpha = transparency
  labs(x = expression(log[2]~fold~change),        # nice maths in labels
       y = expression(-log[10]~italic(p)),
       title = "Volcano plot") +
  theme_bw(base_size = 13)


## 9.2 Mapping vs setting — the distinction that confuses everyone -------
# INSIDE aes(): the property is mapped to a COLUMN of the data,
#               and ggplot creates a legend.
# OUTSIDE aes(): the property is a fixed value for all points.

res_tbl$significant <- res_tbl$padj < 0.1 & abs(res_tbl$log2FC) > 1

ggplot(res_tbl, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(colour = significant), alpha = 0.7) +   # MAPPED -> legend
  scale_colour_manual(values = c("FALSE" = "grey70",
                                 "TRUE"  = "firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  labs(x = expression(log[2]~fold~change), y = expression(-log[10]~italic(p)),
       colour = "FDR<0.1 &\n|log2FC|>1") +
  theme_bw(base_size = 13)

# Compare: colour SET outside aes() — one fixed colour, no legend
ggplot(res_tbl, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(colour = "steelblue", alpha = 0.7) +
  theme_bw()


## 9.3 The geoms you will use this week ----------------------------------

# Histogram — distribution of one continuous variable.
# The p-value histogram is a diagnostic you should ALWAYS look at:
# a flat distribution with a peak near zero is healthy.
ggplot(res_tbl, aes(x = pvalue)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  labs(x = "raw p-value", y = "number of genes") +
  theme_bw()

# Boxplot / violin — a continuous variable split by a category
res_tbl$expr_class <- cut(res_tbl$baseMean,
                          breaks = c(-Inf, 50, 500, Inf),
                          labels = c("low", "medium", "high"))

ggplot(res_tbl, aes(x = expr_class, y = abs(log2FC), fill = expr_class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +   # show the raw data too
  labs(x = "expression level", y = expression(abs(log[2]~FC))) +
  theme_bw() +
  theme(legend.position = "none")

# Bar chart — counts per category
ggplot(res_tbl, aes(x = expr_class, fill = significant)) +
  geom_bar(position = "fill") +          # "fill" = proportions; "stack"; "dodge"
  labs(x = "expression level", y = "proportion of genes") +
  theme_bw()

# Scatter with a trend line — MA plot (mean vs fold change)
ggplot(res_tbl, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(colour = significant), alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = c("grey70", "firebrick")) +
  labs(x = expression(log[10]~mean~expression), y = expression(log[2]~FC)) +
  theme_bw()


## 9.4 Facetting — the same plot, split by a variable --------------------
# One of ggplot's best features. Add one line, get a panel per group.
ggplot(res_tbl, aes(x = log2FC)) +
  geom_histogram(bins = 25, fill = "darkolivegreen4", colour = "white") +
  facet_wrap(~ expr_class, scales = "free_y") +
  labs(x = expression(log[2]~fold~change)) +
  theme_bw()


## 9.5 Saving figures ----------------------------------------------------
# Build the plot as an OBJECT, then save the object. Do this rather than
# relying on the RStudio Plots pane: it is reproducible and lets you
# control size and resolution.

p_volcano <- ggplot(res_tbl, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(colour = significant), alpha = 0.7) +
  scale_colour_manual(values = c("grey70", "firebrick")) +
  labs(x = expression(log[2]~fold~change), y = expression(-log[10]~italic(p))) +
  theme_bw(base_size = 12)

dir.create("figures", showWarnings = FALSE)

ggsave("figures/volcano_demo.png", plot = p_volcano,
       width = 6, height = 4.5, dpi = 300)
# PDF for publication (vector graphics, infinitely scalable):
ggsave("figures/volcano_demo.pdf", plot = p_volcano, width = 6, height = 4.5)

# ALWAYS pass plot = explicitly. A bare ggsave() saves "the last plot",
# which is a silent source of wrong figures in long scripts.


## EXERCISE 9 ---------------------------------------------------------------
# a) Make a histogram of log10(baseMean).
# b) Make a volcano plot where point SIZE reflects log10(baseMean) and
#    colour reflects significance.
# c) Facet your volcano plot by expr_class.
# d) Label the five most significant genes on the volcano plot.
#    Hint: install and use ggrepel::geom_text_repel() with a filtered
#    subset of the data as its own data = argument.
# e) Save your favourite figure as a 300-dpi PNG, 7 x 5 inches.



# ---- 10. Reproducibility: end every script like this ------------------------

# sessionInfo() records R version, platform, and the exact version of every
# loaded package. Paste it into your lab book / methods section. Analyses
# change between package versions; this is how you can explain why.
sessionInfo()

# Consider also: renv::init() to lock package versions per project.


###############################################################################
#
#   WHAT YOU SHOULD BE ABLE TO DO NOW
#   --------------------------------
#   [ ] create and name objects; know their class
#   [ ] build and manipulate vectors; understand vectorisation
#   [ ] compute a log2 fold change from two vectors of counts
#   [ ] load a table, inspect it with str/dim/head, and write it back out
#   [ ] filter rows and select columns, in base R and with dplyr
#   [ ] check that a count matrix and a sample sheet are in the same order
#   [ ] set the reference level of a factor, and say why it matters
#   [ ] write a small function with a default argument
#   [ ] compute per-row and per-column summaries of a matrix
#   [ ] build, modify and save a ggplot figure
#
#   If any box is unticked, that is what to work on in the group session.
#
#   TOMORROW
#   --------
#   The statistics of count data, and your first real differential
#   expression analysis with DESeq2. The filtering, factor-level and
#   column-ordering skills from sections 4, 5 and 8 are the ones you will
#   use immediately.
#
###############################################################################
