# ================================================
# Envrionment setup
# ================================================
getwd()
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
OUT_DIR <- "plots"
dir.create(OUT_DIR, showWarnings = FALSE)

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
  
  sf2 <- read.table(sf2_file, header = TRUE)
  
  p <- ggplot(sf2, aes(x = location / 1e6, y = LR)) +
    geom_point(color = "steelblue", size = 0.8) +
    labs(
      title = population,
      x = "Position (Mb)",
      y = "LR"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[population]] <- p
}

# ================================================
# Combine populations in one column
# ================================================
final_plot <- wrap_plots(plot_list, ncol = 1)

ggsave(
  filename = paste0(OUT_DIR, "/scaffold_", scaffold_id, "_all_populations.png"),
  plot = final_plot,
  width = 8,
  height = 2.5 * length(populations),
  dpi = 300,
  limitsize = FALSE
)
