rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(patchwork)
setwd('~/PhD research/popgen analyses/07. SweepFinder2/polarised')

# ----------------------------------------------------------------
# Plot 1: with invariant sites, polarised (not the entire chromosome - up until ~9.5MB)
IN_DIR <- "sf2_out/A1_AnOudin/"
OUT_DIR <- "plots/A1_AnOudin/"
sf2out1 <- "A1_AnOudin_scaffold_1_invar_wglsfs"
sf2_1 <- read.table(paste0(IN_DIR, sf2out1, ".sf2.out"), header = TRUE)
sf2_1 <- sf2_1[32928:nrow(sf2_1), ]

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
    title = "SweepFinder2 polarised with invariant sites, Chromosome 1, A1_AnOudin",
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

# ----------------------------------------------------------------
# Plot 2: without invariant sites, polarised (the entire chromosome, ~13MB)
IN_DIR <- "sf2_out/A1_AnOudin/"
OUT_DIR <- "plots/A1_AnOudin/"
sf2out2 <- "A1_AnOudin_scaffold_1_wglsfs"
sf2_2 <- read.table(paste0(IN_DIR, sf2out2, ".sf2.out"), header = TRUE)
sf2_2 <- sf2_2[sf2_2$location >= 5850717 & sf2_2$location <= 9419315, ]

# Compute chromosome length from data
chr_length <- max(sf2_2$location, na.rm = TRUE)
width_scale <- chr_length / 1e6   # Mb

# Plot the likelihood ratio (LR) along the chromosome
p2 <- ggplot(sf2_2, aes(x = location, y = LR)) +
  geom_point(color = "steelblue", size = 1) +
#  coord_cartesian(xlim = c(0, chr_length), expand = FALSE) +
  labs(
    title = "SweepFinder2 polarised without invariant sites, Chromosome 1, A1_AnOudin",
    x = "Genomic Position",
    y = "Likelihood Ratio (LR)"
  ) +
  theme_minimal()
p2

# Export the figure
ggsave(paste0(OUT_DIR, sf2out2, "_plot.png"),
       plot = p2,
       width = width_scale,
       height = 4,
       dpi = 300)

# ----------------------------------------------------------------
# Plot 3: unpolarised (not the entire chromosome - up until ~9.5MB)
# the file's path: "C:\Users\u0167113\Documents\PhD research\popgen analyses\07. SweepFinder2\unpolarised\sf2_out\A1_AnOudin_sf2\A1_AnOudin_scaffold_1_wglsfs.sf2.out"
IN_DIR <- "C:/Users/u0167113/Documents/PhD research/popgen analyses/07. SweepFinder2/unpolarised/sf2_out/A1_AnOudin_sf2/"
OUT_DIR <- "plots/A1_AnOudin/"
sf2out3 <- "A1_AnOudin_scaffold_1_wglsfs"
sf2_3 <- read.table(paste0(IN_DIR, sf2out3, ".sf2.out"), header = TRUE)
sf2_3 <- sf2_3[sf2_3$location >= 5850717 & sf2_3$location <= 9419315, ]

# Compute chromosome length from data
chr_length <- max(sf2_3$location, na.rm = TRUE)
width_scale <- chr_length / 1e6   # Mb
max(sf2_3$location, na.rm = TRUE)
min(sf2_3$location, na.rm = TRUE)
# Plot the likelihood ratio (LR) along the chromosome
p3 <- ggplot(sf2_3, aes(x = location, y = LR)) +
  geom_point(color = "steelblue", size = 1) +
  #  coord_cartesian(xlim = c(0, chr_length), expand = FALSE) +
  labs(
    title = "SweepFinder2 unpolarised, Chromosome 1, A1_AnOudin",
    x = "Genomic Position",
    y = "Likelihood Ratio (LR)"
  ) +
  theme_minimal()
p3

# Export the figure
ggsave(paste0(OUT_DIR, sf2out3, "unpolarised_plot.png"),
       plot = p3,
       width = width_scale,
       height = 4,
       dpi = 300)

# ----------------------------------------------------------------
# Combine plots

combined <- p3 / p1 / p2  # '/' stacks, '|' would place side-by-side
combined

ggsave(
  file.path(OUT_DIR, "comparison.png"),
  plot = combined,
  width = width_scale,
  height = 12,
  dpi = 300
)

