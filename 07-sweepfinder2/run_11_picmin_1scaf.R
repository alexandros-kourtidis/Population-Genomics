# ==========================================
# PicMin
# https://github.com/TBooker/PicMin/tree/main
# ==========================================

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
# 1) User settings
# ------------------------------------------
# - Specify: input/output files, population column names,
#    window size, target scaffold, and PicMin simulation parameters.
# - The scaffold is set as an argument from a SLURM array.
#    Example bash SLURM array command:
#     SC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffolds.txt)
#     Rscript run_picmin_1scaff.R $SC
#    (you can also manually define it yourself and comment out the args variable assignment)
# - For missing_levels: nmax=npopcols. Avoid n = 2 (low power)
# ==========================================
IN_DIR		<- "input"
OUT_DIR		<- "output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pop_cols	<- c("C1_BlfN", "C2_MO", "C3_Ter1", "C4_BW_48630", "C5_BW_36962")
missing_levels 	<- c(3, 4, 5)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("ERROR: No scaffold argument supplied. Usage: Rscript run_picmin_1scaff.R <scaffold>")
}
SCAFFOLD <- args[1]

WIN_SIZE	<- 10000L
numReps		<- 100000	# increase to 1e5 for final accuracy
nsims		<- 40000	# increase to 4e4 for final accuracy

infile		<- file.path(IN_DIR, "PicMin_input_meanLR_10kb_windows_woContam.txt")
outfile 	<- file.path(OUT_DIR, paste0("popsC1C2C3C4C5_", SCAFFOLD, ".tsv"))

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
# ensuring row names match PicMin’s required locus ide
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
