### ============================================
### Load libraries
### ============================================
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)

rm(list = ls())

### ============================================
### Set your PIXY output directory
### ============================================
pixy_dir <- "."   # <-- change if needed
out <- "plots"

### ============================================
### Function to read all PIXY files for a stat
### ============================================
read_pixy_stat <- function(stat, subfolder = "pixy_stats") {
  stat_dir <- file.path(pixy_dir, subfolder)
  files <- list.files(
    stat_dir,
    pattern = paste0(stat, ".txt$"),
    full.names = TRUE
  )
  dt <- data.table::rbindlist(lapply(files, data.table::fread))
  dt$stat <- stat
  return(dt)
}

### ============================================
### Load fst stats
### ============================================
fst_df    <- read_pixy_stat("fst", "pixy_stats")

# ============================================
# Select the B2vsC1 comparison
# ============================================
fst_pair <- fst_df %>%
  filter(
    (pop1 == "B2_OM2" & pop2 == "C1_BlfN") |
      (pop1 == "C1_BlfN" & pop2 == "B2_OM2")
  )

# ============================================
# Set the genome axis
# ============================================
chr_lengths <- fst_pair %>%
  group_by(chromosome) %>%
  summarise(chr_len = max(window_pos_2), .groups = "drop") %>%
  mutate(
    chr_num = as.numeric(stringr::str_extract(chromosome, "\\d+"))
  ) %>%
  arrange(chr_num)

chr_lengths <- chr_lengths %>%
  mutate(
    chr_start = lag(cumsum(chr_len), default = 0)
  )

fst_pair <- fst_pair %>%
  left_join(chr_lengths, by = "chromosome") %>%
  mutate(
    window_mid = (window_pos_1 + window_pos_2) / 2,
    genome_pos = chr_start + window_mid
  )

chr_labels <- chr_lengths %>%
  mutate(
    label_pos = chr_start + chr_len / 2
  )

# ============================================
# Set the colours
# ============================================
chr_lengths <- chr_lengths %>%
  mutate(
    chr_index = row_number(),
    chr_shade = if_else(chr_index %% 2 == 1, "dark", "light")
  )
fst_pair <- fst_pair %>%
  left_join(
    chr_lengths %>% select(chromosome, chr_shade),
    by = "chromosome"
  )
# ============================================
# Plot Fst accross the genomehttp://127.0.0.1:8497/graphics/plot_zoom_png?width=1904&height=1090
# ============================================

ggplot(fst_pair, aes(x = genome_pos, y = avg_wc_fst, color = chr_shade)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_x_continuous(
    breaks = chr_labels$label_pos,
    labels = chr_labels$chromosome,
    expand = c(0.01, 0)
  ) +
  scale_color_manual(
    values =c(dark = "grey30", light = "grey70")
  ) +
  labs(
    title = "Genome-wide FST: B2_OM2 vs C1_BlfN",
    x = "Genomic position",
    y = "Weir & Cockerham FST"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)
  )


# ============================================
# Plot Fst scans per chromosome
# ============================================

p_pair <- ggplot(fst_pair, aes(x = window_mid, y = avg_wc_fst)) +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 4) +
  theme_bw() +
  labs(
    x = "Position within chromosome",
    y = "FST"
  )

ggsave(file.path(out, "fst_scan_B2_C1.png"), p_pair, dpi = 300)

# ============================================
# Select the B2vsC1 comparison
# ============================================
# ============================================
# Select the B2vsC1 comparison
# ============================================