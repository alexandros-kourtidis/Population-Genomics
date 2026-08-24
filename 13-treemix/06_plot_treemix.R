## ============================================================
## Compare TreeMix population tree topologies across m values
## ============================================================
## Assumes TreeMix was run per-m as, e.g.:
##   treemix -i data.treemix.gz -m 4 -o treemix_mruns/treemix.m4
## i.e. output files look like: treemix_mruns/treemix.m4.treeout.gz, .cov.gz, etc.
## Adjust `base_prefix` below if your naming differs.

# --- Paths ---
treemix_dir <- "treemix_mruns"   # folder with TreeMix output files
plots_dir   <- "plots"           # folder to save plots
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

base_prefix <- "treemix"         # prefix used before ".mX" in your output filenames

# m values to compare
m_values <- c(0, 1, 2, 4, 6)

# --- Get TreeMix's plotting helper functions (downloads once if missing) ---
source("/vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/envs/treemix/bin/plotting_funcs.R")

# Optional: file listing population order, one pop per line (used for residual plots)
# Set to NULL if you don't have one / don't want residual plots.
poporder_file <- "pop_list.txt"

# Read the UTF-16 file properly, then save a clean UTF-8 version
con <- file("pop_list.txt", encoding = "UTF-16LE")
pops <- readLines(con)
close(con)

pops <- pops[pops != ""]   # drop any blank lines
writeLines(pops, "poporder_utf8.txt")

# check it looks right
pops

# ------------------------------------------------
# Test plot for m0
# ------------------------------------------------
# the stem of the files for this m
stem <- file.path(treemix_dir, paste0(base_prefix, ".m0.rep1"))

# plot tree
p<- plot_tree(stem)
p
title(main = "TreeMix tree, m = 0")

# residuals
plot_resid(stem, poporder_file)
title(main = "Residuals, m = 0")

# plot it
png("test_tree_m0.png", width = 1600, height = 1200, res = 200)
plot_tree(stem)
title(main = "TreeMix tree, m = 0")
dev.off()

# ------------------------------------------------
# Per-m plots: tree + residuals, saved individually
# ------------------------------------------------
for (m in m_values) {
  stem <- file.path(treemix_dir, paste0(base_prefix, ".m", m, ".rep1"))
  
  if (!file.exists(paste0(stem, ".treeout.gz"))) {
    warning("Skipping m=", m, ": expected file not found: ", stem, ".treeout.gz")
    next
  }
  
  message("Plotting m = ", m)
  
  # Tree topology
  png(file.path(plots_dir, paste0("tree_m", m, ".png")),
      width = 1600, height = 1200, res = 200)
  plot_tree(stem)
  title(main = paste0("TreeMix tree, m = ", m))
  dev.off()
  
  # Residuals (how well the tree + m migration edges fit the covariance)
  if (!is.null(poporder_file)) {
    png(file.path(plots_dir, paste0("resid_m", m, ".png")),
        width = 1600, height = 1200, res = 200)
    plot_resid(stem, poporder_file)
    title(main = paste0("Residuals, m = ", m))
    dev.off()
  }
}

# --- Combined grid figure: all trees side-by-side for quick topology comparison ---
valid_m <- Filter(function(m) file.exists(file.path(treemix_dir, paste0(base_prefix, ".m", m, ".rep1", ".treeout.gz"))), m_values)

if (length(valid_m) > 0) {
  n <- length(valid_m)
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  
  png(file.path(plots_dir, "trees_comparison_grid.png"),
      width = 600 * ncol, height = 500 * nrow, res = 150)
  par(mfrow = c(nrow, ncol))
  for (m in valid_m) {
    stem <- file.path(treemix_dir, paste0(base_prefix, ".m", m))
    plot_tree(stem)
    title(main = paste0("m = ", m))
  }
  dev.off()
  message("Combined comparison grid saved: ", file.path(plots_dir, "trees_comparison_grid.png"))
}

message("Done. Individual and combined plots saved in: ", normalizePath(plots_dir))