rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(patchwork)

# Define populations
populations <- c("A1_AnOudin", "A2_BuSN", "A3_Mech", "A4_PL15_YEL", "A5_GenM",
                 "B1_BKN1", "B2_OM2", "B3_ZW", "B4_OHZ", "B5_DA2",
                 "C1_BlfN", "C2_MO", "C3_Ter1", "C4_BW_48630", "C5_BW_36962",
                 "D1_CBOO6", "D2_LRV", "D3_BKLE5", "D4_BW_62256", "D5_BW_22050")

# Define background colors by prefix
bg_colors <- list(
  "A" = "#e6f9e6",  # light green
  "B" = "#e6f0ff",  # light blue
  "C" = "#ffffe6",  # light yellow
  "D" = "#ffe6e6"   # light red
)

# Create output directory
OUT_DIR <- "plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Initialize list to store combined plots per population
all_population_plots <- list()

for (population in populations) {
  IN_DIR <- paste0("sf2_out/",population)
  
  scaffold_files <- file.path(IN_DIR, paste0(population, "_scaffold_", 1:10, "_wglsfs.sf2.out"))
  plot_list <- list()
  
  # Determine global y-axis limits for this population
  y_min <- Inf
  y_max <- -Inf
  for (i in 1:10) {
    sf2gl <- read.table(scaffold_files[i], header = TRUE)
    y_min <- min(y_min, min(sf2gl$LR, na.rm = TRUE))
    y_max <- max(y_max, max(sf2gl$LR, na.rm = TRUE))
  }
  
  # Get background color based on population prefix
  prefix <- substr(population, 1, 1)
  bg_color <- bg_colors[[prefix]]
  
  # Generate plots for each scaffold
  for (i in 1:10) {
    sf2gl <- read.table(scaffold_files[i], header = TRUE)
    p <- ggplot(sf2gl, aes(x = location, y = LR)) +
      geom_point(color = "steelblue", size = 1) +
      labs(title = NULL, x = NULL, y = NULL) +
      ylim(y_min, y_max) +
      theme_minimal(base_size = 8) +
      theme(
        plot.background = element_rect(fill = bg_color, color = NA),
        panel.background = element_rect(fill = bg_color, color = NA)
      )
    plot_list[[i]] <- p
  }
  
  # Combine scaffold plots horizontally for this population
  combined_population_plot <- wrap_plots(plot_list, nrow = 1)
  all_population_plots[[population]] <- combined_population_plot
}

# Combine all population plots vertically
final_plot <- wrap_plots(all_population_plots, ncol = 2)

# Save the final combined plot
ggsave(filename = paste0(OUT_DIR, "/all_populations_2columns.png"), plot = final_plot,
       width = 30, height = 4 * length(populations), dpi = 300, limitsize = FALSE)

# Save each population's plot in a separate file
for (i in 1:10) {
  ggsave(
    filename = paste0(OUT_DIR, "/", population, "_scaffold_", i, ".png"),
    plot = plot_list[[i]],
    width = 8, height = 4, dpi = 300
  )
}


# Save each population's plot
names(all_population_plots) <- populations

for (pop in names(all_population_plots)) {
  ggsave(
    filename = file.path(OUT_DIR,pop, paste0(pop, "_combined.png")),
    plot = all_population_plots[[pop]],
    width = 30, height = 4, dpi = 300, limitsize = FALSE
  )
}

# To test-save
# ggsave(filename = paste0(OUT_DIR, "/all_populations_combined_sf2_nolabels_colours_test.png"), plot = final_plot, width = 10, height = 1 * length(populations), dpi = 50, limitsize = FALSE)
