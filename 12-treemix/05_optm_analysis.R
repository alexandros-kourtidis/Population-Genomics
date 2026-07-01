# OptM analysis for TreeMix outputs
# Expects a single flat folder containing all .llik / .cov.gz / .modelcov.gz

rm(list=ls())

library(OptM)

indir <- "treemix_mruns"
outdir <- "OptM_results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


folder <- "treemix_mruns"

# Default method
test.optM <- optM(indir)
plot_optM(test.optM, method = "Evanno")

# Linear/change-point methods
test.linear <- optM(indir, method = "linear")
plot_optM(test.linear, method = "linear")

# SiZer method
test.sizer <- optM(indir, method = "SiZer")
plot_optM(test.sizer, method = "SiZer")
