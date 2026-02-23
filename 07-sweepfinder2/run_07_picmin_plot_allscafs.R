# ==========================================
# Plotting PicMin results
# ------------------------------------------
# Create a Manhattan-style scatter plot showing signal strength (-log10 q)
# along the scaffold, colored by the estimated number of contributing lineages.
# ! You can FILTER to plot only for the desire number of lineages !
# ==========================================
rm(list = ls())
library(data.table)
library(ggplot2)
library(cowplot)

setwd("C:/Users/u0167113/Documents/PhD research/popgen analyses/07. SweepFinder2/polarised/picmin")

# ================================================
# User settings
# ================================================
IN_DIR  <- "results/pops_agriculture"
OUT_DIR   <- "plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

input_pops <- "pops_agriculture_"
scaffold_files <- file.path(IN_DIR, paste0(input_pops, "scaffold_", 1:10, ".tsv"))
output <- "pops_agriculture"
  
# Load & combine
df_list <- lapply(1:10, function(i){
  f <- scaffold_files[i]
  d <- fread(f)
  d[, scaffold := i]          # numeric scaffold ID for sorting
  d[, scaffold_lab := paste0("Scaffold ", i)]
  return(d)
})
df_all_alln <- rbindlist(df_list)

# ================================================
# FILTER : keep only lineages more than []
nlin <- 7
df_all <- df_all_alln[n_est > nlin]
# ================================================

# Manually force scaffold order 1 → 10
df_all[, scaffold_lab := factor(scaffold_lab, levels = paste0("Scaffold ", 1:10))]

# Color palette - rainbow gradient
col_pal <- c(
  "#0000FF", #"#0033FF",
  "#0066FF", #"#0099FF",
  "#00CCFF", #"#00FFEE",
  "#00FFBB", #"#00FF88",
  "#00FF55", #"#33FF22",
  "#66FF00", #"#99FF00",
  "#CCFF00", #"#FFFF00",
  "#FFCC00", #"#FF9900",
  "#FF6600", #"#FF3300",
  "#FF1900", #"#FF0000"  
)

sig_line <- -log10(0.05)

# Alternating facet background colors
bg_df <- data.table(
  scaffold_lab = factor(paste0("Scaffold ", 1:10), levels = paste0("Scaffold ", 1:10)),
  bg_col       = rep(c("white", "gray95"), length.out = 10)
)

# Determine shared y‑axis range
YMAX <- max(-log10(df_all$pooled_q), na.rm = TRUE)

# ================================================
#  Make the plot
# ================================================
p <- ggplot(
  df_all,
  aes(
    x = start / 1e6,
    y = -log10(pooled_q),
    fill = factor(n_est)
  )
) +
  # Panel background rectangles
  geom_rect(
    data = bg_df,
    inherit.aes = FALSE,
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf,
    ),
    fill = bg_df$bg_col,
    alpha = 0.4
  ) +
  # Actual SNP points
  geom_point(shape = 21, size = 2, alpha = 0.3) +
  geom_hline(yintercept = sig_line, lty = 2) +
  
  # --- Important: fill scale only for n_est, not for backgrounds!
  scale_fill_manual(
    name   = "Number of Lineages",
    values = col_pal,
    guide  = guide_legend(order = 1)
  ) +
  
  scale_y_continuous(
    expression(-log[10] * "(q-value)"),
    limits = c(0, YMAX)
  ) +
  scale_x_continuous("Position (Mbp)",
                     breaks = seq(0, 20, by = 3)) +
  
  facet_wrap(~ scaffold_lab, 
             nrow = 1, 
             strip.position = "bottom",
             scales = "free_x" ) +
  
  theme_half_open() +
  background_grid() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "bottom",
    legend.margin = margin(t = -5),
    legend.title = element_text(size = 10),
    
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    
    panel.spacing = unit(0.2, "lines"),
    axis.text.x = element_text(size = 6)
  )

# Add the title
p <- p +
  ggtitle(input_pops) +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))


print(p)

# ================================================
# Save PicMin results into high-quality PNG
# ================================================
png_file <- file.path(
  OUT_DIR,
  paste0(output, "_above", nlin, "lines_allscafs_plot.png")
)

ggsave(
  filename = png_file,
  plot     = p,
  width    = 36,
  height   = 6,
  dpi      = 300,
  bg       = "white",
  device   = grDevices::png,
  type     = "cairo"
)

message("Saved: ", png_file)
