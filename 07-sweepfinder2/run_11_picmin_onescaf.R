# ==========================================================
# PicMin
# https://github.com/TBooker/PicMin/tree/main
# ==========================================================
# PicMin R Script – SLURM-Controlled Parameter Interface
# ----------------------------------------------------------
# This script is designed to be run from a SLURM array job. All required input parameters are supplied by the SLURM
# submission script using command-line arguments. Nothing is hard-coded inside this R script.
#
# The script expects EXACTLY 10 arguments in the following order (matching the SLURM script):
#
#   args[1] = SCAFFOLD
#       Name of the genomic scaffold to process. The SLURM
#       array determines which scaffold is passed.
#
#   args[2] = pop_cols
#       Comma-separated list of population column names
#       present in the input file (no spaces). This vector
#       is split in R and defines the population-specific
#       p-values used by PicMin.
#
#   args[3] = missing_levels
#       Comma-separated integers specifying how many
#       populations must be present (non-missing) for each
#       PicMin run. For each n in this list, a separate
#       null model is generated and applied.
#
#   args[4] = IN_DIR
#       Directory where the input file is located.
#
#   args[5] = OUT_DIR
#       Directory where output files will be written.
#       Created automatically if not present.
#
#   args[6] = infile
#       Full path to the PicMin input file containing all
#       scaffolds. The script subsets it based on SCAFFOLD.
#
#   args[7] = OUTPREFIX
#       Prefix used when constructing the output filename.
#       The final output file will be:
#            <OUT_DIR>/<OUTPREFIX><SCAFFOLD>.tsv
#
#   args[8] = WIN_SIZE
#       Window size used by the input dataset (e.g. 10000).
#       This does not change PicMin internals but is kept
#       for clarity and future flexibility.
#
#   args[9] = numReps
#       Number of replicates used by PicMin for estimating
#       the p-value for each window. Higher = more accurate.
#
#  args[10] = nsims
#       Number of null simulations used to generate order-
#       statistic correlation matrices. Larger values improve
#       stability but increase memory usage.
#
# ----------------------------------------------------------
# GENERAL NOTES:
# - All argument parsing is strict; the script halts if the
#   expected number of parameters is not supplied.
# - The script constructs its output filename internally
#   using OUTPREFIX and SCAFFOLD.
# - No assumptions are made about job arrays; they are
#   handled entirely by the SLURM script.
# - Null correlation matrices are generated per missing-
#   level category, ensuring correct PicMin inference.
# ==========================================================

rm(list = ls())

# ==========================================
# 0) Environment setup
# ==========================================
library(data.table)
library(dplyr)
library(PicMin)
library(poolr)

# install.packages("remotes")
# remotes::install_github("TBooker/PicMin")
# install.packages(c("tidyverse","poolr","cowplot"))

setwd("/lustre1/scratch/363/vsc36396/picmin")

# ==========================================
# 1) USER SETTINGS
# ==========================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
    stop("Not enough arguments supplied. Expected 10 arguments: SCAFFOLD, POPS, MISS, IN_DIR, OUT_DIR, INFILE, OUTPREFIX, WINDOW, NUMREPS, NSIMS")
}

SCAFFOLD       <- args[1]
pop_cols       <- unlist(strsplit(args[2], ","))
missing_levels <- as.integer(strsplit(args[3], ",")[[1]])
IN_DIR         <- args[4]
OUT_DIR        <- args[5]
infile         <- args[6]
OUTPREFIX      <- args[7]
WIN_SIZE       <- as.integer(args[8])
numReps        <- as.integer(args[9])
nsims          <- as.integer(args[10])

outfile <- file.path(OUT_DIR, paste0(OUTPREFIX, SCAFFOLD, ".tsv"))
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# START MESSAGE
message(">>> PicMin analysis started at: ", Sys.time())
message(">>> Running on scaffold: ", SCAFFOLD)
message(">>> Populations included: ", paste(pop_cols, collapse = ", "))
message(">>> Using numReps = ", numReps, " and nsims = ", nsims)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ==========================================
# 2) Read file & subset to scaffold
# ==========================================
x <- fread(infile)
x <- x[scaffold == SCAFFOLD]

# ==========================================
# 3) Compute empirical p-values
# ==========================================
stats <- x[, ..pop_cols]

pvals <- as.data.table(
  lapply(
    stats,
    function(col) PicMin:::EmpiricalPs(col, large_i_small_p = TRUE)
  )
)

setnames(pvals, pop_cols)
pvals[, window_id := x$window_id]

# ==========================================
# 4) Build p-value matrix for PicMin
# ------------------------------------------
# Convert p-values to a matrix with windows as rows and lineages as columns,
# ensuring row names match PicMinâ€™s required locus ide
# ==========================================
pmat <- as.matrix(pvals[, ..pop_cols])
rownames(pmat) <- pvals$window_id
nLins_total <- ncol(pmat)

# ==========================================
# 5) Build null matrices & run PicMin
# ------------------------------------------
# For each missing-data level (3/4/5 lineages present):
# simulate null p-values, compute null correlations,
# then run PicMin to obtain p- and q-values per window.
# ==========================================
results_list <- list()

for (n in missing_levels) {
  
  message("=== Running PicMin for n = ", n, " populations present ===")
  
  # Null correlation of ordered p-values
  null_raw     <- t(replicate(nsims, PicMin:::GenerateNullData(1.0, n, 0.5, 3, 10000)))
  null_ordered <- t(apply(null_raw, 1, PicMin:::orderStatsPValues))
  null_cor     <- cor(null_ordered)
  
  # Windows with exactly n non-missing populations
  keep <- rowSums(is.na(pmat)) == (nLins_total - n)
  if (!any(keep)) next
  
  pmat_n <- pmat[keep, , drop = FALSE]
  
  # Run PicMin
  res_p <- res_n <- numeric(nrow(pmat_n))
  
  for (i in seq_len(nrow(pmat_n))) {
    pvec <- na.omit(pmat_n[i, ])
    tr   <- PicMin:::PicMin(pvec, null_cor, numReps = numReps)
    
    res_p[i] <- tr$p
    res_n[i] <- tr$config_est
  }
  
  df_out <- data.frame(
    numLin = n,
    p      = res_p,
    q      = p.adjust(res_p, "fdr"),
    n_est  = res_n,
    locus  = rownames(pmat_n),
    stringsAsFactors = FALSE
  )
  
  results_list[[paste0("n", n)]] <- df_out
}

# ==========================================
# 6) Merge results with coordinates
# ------------------------------------------
# Add genomic coordinates (scaffold, start, end) to each window,
# enabling positional interpretation and plotting.
# ==========================================
picmin_results <- dplyr::bind_rows(results_list)
picmin_results$pooled_q <- p.adjust(picmin_results$p, "fdr")  
meta <- x[, .(window_id, scaffold, start, end)]

picmin_results <- merge(
  as.data.table(picmin_results),
  as.data.table(meta),
  by.x = "locus",
  by.y = "window_id",
  all.x = TRUE
)

# ==========================================
# 7) Save output
# ------------------------------------------
# Write the full PicMin results table to file:
# p/q-values, estimated lineage counts, and coordinates.
# ==========================================
fwrite(picmin_results, file = outfile, sep = "\t")
message("head(picmin_results):")
print(head(picmin_results))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# END MESSAGE
message(">>> PicMin analysis COMPLETED at: ", Sys.time())
message(">>> Output written to: ", outfile)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
