# OptM analysis for TreeMix outputs
# Expects a single flat folder containing all .llik / .cov.gz / .modelcov.gz

rm(list=ls())

library(OptM)

indir <- "treemix_mruns"
outdir <- "OptM_results/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Default method
test.optM <- optM(indir)
plot_optM(test.optM, method = "Evanno")


# Linear/change-point methods
test.linear <- optM(indir, method = "linear")
plot_optM(test.linear, method = "linear")

png(filenam=paste0(outdir, "OptM_linear.png"), 
    width = 2000,
    height = 1700,
    res = 300)
plot_optM(test.linear, method = "linear")
dev.off()

# SiZer method
test.sizer <- optM(indir, method = "SiZer")
plot_optM(test.sizer, method = "SiZer")

png(filenam=paste0(outdir, "OptM_SiZer.png"), 
    width = 2000,
    height = 1700,
    res = 300)
plot_optM(test.sizer, method = "SiZer")
dev.off()

