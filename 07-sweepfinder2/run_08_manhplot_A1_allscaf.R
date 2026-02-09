rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(patchwork)  # install.packages("patchwork") if needed

# ---------------------------------------------------------------
# Define population and directories
setwd('~/PhD research/popgen analyses/07. SweepFinder2/polarised')

population <- "A1_AnOudin"

OUT_DIR <- paste0("plots/", population)
IN_DIR <- paste0("sf2_out/", population)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

scaffold_files <- file.path(IN_DIR, paste0(population, "_scaffold_", 1:10, "_wglsfs.sf2.out"))
plot_list <- list()

# ---------------------------------------------------------------
# Make the plots

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
    geom_point(color = "steelblue", size = 0.1) +
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
ggsave(filename = file.path(OUT_DIR, paste0(population, "_AllScaffolds.sf2.png")), combined_plot, width = 30, height = 4, dpi = 300)
