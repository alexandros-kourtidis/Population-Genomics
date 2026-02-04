# Load necessary libraries
#library(ggplot2)

# Read the SweepFinder2 output file
#sf2 <- read.table("A1_AnOudin_scaffold_10.sf2.out.txt", header = TRUE)
#
# Plot the likelihood ratio (LR) along the chromosome
#ggplot(sf2, aes(x = location, y = LR)) +
#  geom_point(color = "steelblue", size = 1) +
#  labs(
#    title = "SweepFinder2 on Chromosome 10, A1_AnOudin",
#    x = "Genomic Position",
#    y = "Likelihood Ratio (LR)"
#  ) +
#  theme_minimal()

# Read the SweepFinder2 WITH GLOBAL SFS output file
#sf2gl <- read.table("A1_AnOudin_scaffold_10_wglsfs.sf2.out.txt", header = TRUE)

# Plot the likelihood ratio (LR) along the chromosome
#ggplot(sf2gl, aes(x = location, y = LR)) +
#  geom_point(color = "steelblue", size = 1) +
#  labs(
#    title = "SweepFinder2 with global sfs, Chromosome 10, A1_AnOudin",
#    x = "Genomic Position",
#    y = "Likelihood Ratio (LR)"
#  ) +
#  theme_minimal()


rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(patchwork)  # install.packages("patchwork") if needed

# Define population and directories
population <- "A1_AnOudin"

OUT_DIR <- paste0(population, "_plots")
IN_DIR <- paste0(population, "_sf2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# List of scaffold file names and empty list for the plots
scaffold_files <- file.path(IN_DIR, paste0(population, "_scaffold_", 1:10, "_wglsfs.sf2.out"))
plot_list <- list()


# Determine global y-axis limits
y_min <- Inf
y_max <- -Inf

for (i in 1:10) {
  sf2gl <- read.table(scaffold_files[i], header = TRUE)
  y_min <- min(y_min, min(sf2gl$LR, na.rm = TRUE))
  y_max <- max(y_max, max(sf2gl$LR, na.rm = TRUE))
}

# Loop through each scaffold file
for (i in 1:10) {
  sf2gl <- read.table(scaffold_files[i], header = TRUE)
  
  p <- ggplot(sf2gl, aes(x = location, y = LR)) +
    geom_point(color = "steelblue", size = 1) +
    labs(
      title = paste("SF2 ", population, ", Chromosome", i),
      x = "Genomic Position",
      y = NULL
    ) +
    ylim(y_min, y_max) +
    theme_minimal(base_size = 8)
  
  # Save individual plot
  ggsave(filename = file.path(OUT_DIR, paste0(population, "_Scaffold", i, ".sf2.png")), plot = p, width = 8, height = 4, dpi = 300)
  
  # Store plot in list
  plot_list[[i]] <- p
}

# Combine all plots horizontally
combined_plot <- wrap_plots(plot_list, nrow = 1) &
  theme(axis.title.y = element_text("Likelihood Ratio (LR)"))

# Save the combined plot
ggsave(filename = file.path(OUT_DIR, paste0("combined_sweepfinder2_plots_", population, "sf2.png")), combined_plot, width = 30, height = 4, dpi = 300)
