# ================================================
# Environment setup
# ================================================
rm(list = ls())

library(ggplot2)
library(patchwork)

setwd("~/PhD research/popgen analyses/07. SweepFinder2/polarised")

# ================================================
# User settings
# ================================================
populations <- c(
  "A1_AnOudin", "A2_BuSN", "A3_Mech", "A4_PL15_YEL", "A5_GenM",
  "B1_BKN1", "B2_OM2", "B3_ZW", "B4_OHZ", "B5_DA2",
  "C1_BlfN", "C2_MO", "C3_Ter1", "C4_BW_48630", "C5_BW_36962",
  "D1_CBOO6", "D2_LRV", "D3_BKLE5", "D4_BW_62256", "D5_BW_22050"
)

scaffold_id <- 5
region_mb <- c(0.9, 1.4)

OUT_DIR <- "plots/region_0.9_1.4Mb_noglobalyaxis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_list <- list()

# ================================================
# Loop over populations
# ================================================
for (population in populations) {
  
  message("Plotting: ", population)
  
  sf2_file <- file.path(
    "sf2_out",
    population,
    paste0(population, "_scaffold_", scaffold_id, "_wglsfs.sf2.out")
  )
  
  # Read full scaffold
  sf2 <- read.table(sf2_file, header = TRUE)
  
  # ---- Y-axis based on ENTIRE scaffold ----
  y_min <- min(sf2$LR, na.rm = TRUE)
  y_max <- max(sf2$LR, na.rm = TRUE)
  
  padding <- 0.05 * (y_max - y_min)
  y_limits <- c(y_min - padding, y_max + padding)
  
  # ---- Subset zoomed region ----
  sf2_sub <- subset(
    sf2,
    location >= region_mb[1] * 1e6 &
      location <= region_mb[2] * 1e6
  )
  
  p <- ggplot(sf2_sub, aes(x = location / 1e6, y = LR)) +
    geom_point(color = "steelblue", size = 0.8) +
    coord_cartesian(
      xlim = region_mb,
      ylim = y_limits
    ) +
    labs(
      title = population,
      x = "Position (Mb)",
      y = "LR"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  
  # Save individual plot
  ggsave(
    filename = paste0(
      OUT_DIR, "/",
      population,
      "_scaffold_", scaffold_id,
      "_", region_mb[1], "_", region_mb[2], "Mb.png"
    ),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300,
    bg = "white"
  )
  
  plot_list[[population]] <- p
}

# ================================================
# Combine all populations vertically
# ================================================
final_plot <- wrap_plots(plot_list, ncol = 1)

ggsave(
  filename = paste0(
    OUT_DIR, "/",
    "scaffold_", scaffold_id,
    "_", region_mb[1], "_", region_mb[2],
    "Mb_all_populations_vertical.png"
  ),
  plot = final_plot,
  width = 7,
  height = 2.3 * length(populations),
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)
