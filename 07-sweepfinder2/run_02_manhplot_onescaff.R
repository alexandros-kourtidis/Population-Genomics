rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(patchwork)
setwd('~/PhD research/popgen analyses/07. SweepFinder2/polarised')

# ----------------------------------------------------------------
# Plot 1: AnOudin Scaffold 1 WITH global sfs
IN_DIR <- "sf2_out/A1_AnOudin/"
OUT_DIR <- "plots/A1_AnOudin/"
sf2out1 <- "A1_AnOudin_scaffold_1_wglsfs"
sf2_1 <- read.table(paste0(IN_DIR, sf2out1, ".sf2.out"), header = TRUE)

# Compute chromosome length from data
chr_length <- max(sf2_1$location, na.rm = TRUE)
width_scale <- chr_length / 1e6   # Mb
max(sf2_1$location, na.rm = TRUE)
min(sf2_1$location, na.rm = TRUE)
# Plot the likelihood ratio (LR) along the chromosome
p1 <- ggplot(sf2_1, aes(x = location, y = LR)) +
  geom_point(color = "steelblue", size = 1) +
  #  coord_cartesian(xlim = c(0, chr_length), expand = FALSE) +
  labs(
    title = "SweepFinder2 polarised, Chromosome 1, A1_AnOudin",
    x = "Genomic Position",
    y = "Likelihood Ratio (LR)"
  ) +
  theme_minimal()
p1

# Export the figure
ggsave(paste0(OUT_DIR, sf2out1, "_plot.png"),
       plot = p1,
       width = width_scale,
       height = 4,
       dpi = 300)
