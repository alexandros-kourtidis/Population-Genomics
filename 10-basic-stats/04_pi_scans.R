### ============================================
### Load libraries
### ============================================
library(tidyverse)
library(data.table)

rm(list = ls())

### ============================================
### Set directories
### ============================================
pixy_dir <- "."
out <- "plots"

### ============================================
### Function to read PIXY stats
### ============================================
read_pixy_stat <- function(stat, subfolder = "pixy_stats") {
  stat_dir <- file.path(pixy_dir, subfolder)
  files <- list.files(
    stat_dir,
    pattern = paste0(stat, ".txt$"),
    full.names = TRUE
  )
  dt <- rbindlist(lapply(files, fread))
  dt$stat <- stat
  return(dt)
}

### ============================================
### Load pi stats
### ============================================
pi_df <- read_pixy_stat("pi", "pixy_stats")

# Optional: remove NA windows
pi_df <- pi_df %>% filter(!is.na(avg_pi))

# Optional: filter low-SNP windows
pi_df <- pi_df %>% filter(no_snps >= 10)

# ============================================
# Order scaffolds & set x axis
# ============================================

chr_lengths <- pi_df %>%
  group_by(chromosome) %>%
  summarise(chr_len = max(window_pos_2), .groups = "drop") %>%
  mutate(
    chr_num = as.numeric(stringr::str_extract(chromosome, "\\d+"))
  ) %>%
  arrange(chr_num) %>%
  mutate(
    chr_start = lag(cumsum(chr_len), default = 0),
    chr_index = row_number(),
    chr_shade = if_else(chr_index %% 2 == 1, "dark", "light")
  )

pi_df <- pi_df %>%
  left_join(
    chr_lengths %>% select(chromosome, chr_start, chr_shade),
    by = "chromosome"
  ) %>%
  mutate(
    window_mid = (window_pos_1 + window_pos_2) / 2,
    genome_pos = chr_start + window_mid
  )

chr_labels <- chr_lengths %>%
  mutate(label_pos = chr_start + chr_len / 2)

# ============================================
# Plot all pops
# ============================================

p_all <- ggplot(pi_df, aes(x = genome_pos, y = avg_pi, color = chr_shade)) +
  geom_point(size = 0.35, alpha = 0.6) +
  scale_color_manual(
    values = c(dark = "grey30", light = "grey70")
  ) +
  scale_x_continuous(
    breaks = chr_labels$label_pos,
    labels = chr_labels$chromosome,
    expand = c(0.01, 0)
  ) +
  labs(
    x = "Genomic position",
    y = expression(pi),
    title = "Genome-wide nucleotide diversity (π)"
  ) +
  facet_wrap(~ pop, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 7),
    strip.text = element_text(face = "bold")
  )
p_all

# Save them all in a single png
ggsave(
  filename = file.path(out,"pi_genomewide_all_populations.png"),
  plot = p_all,
  width = 14,
  height = 10,
  units = "in",
  dpi = 300
)
