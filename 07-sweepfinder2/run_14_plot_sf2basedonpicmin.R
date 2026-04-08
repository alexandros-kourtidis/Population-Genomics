# ==========================================
# Plotting SF2 data based on PicMin results
# ------------------------------------------
# It takes your two inputs (the 500 bp SweepFinder2 output + the 10 kb PicMin q‑values)
# and produces a merged table where each 500 bp window receives a label (“selected” if 
# it falls inside any PicMin 10 kb region with q < 0.05; otherwise “neutral”).
#
# !!! IMPORTANT: now it takes all number of lineages from picmin, we should probably control that
# ==========================================
rm(list = ls())
library(dplyr)
library(readr)
library(ggplot2)

setwd("C:/Users/u0167113/Documents/PhD research/popgen analyses/07. SweepFinder2/polarised")

# =============================
# 0. Directories
# =============================

IN_DIR_SF2 <- "sf2_out/C1_BlfN"
IN_DIR_Pic <- "picmin/results/pops_C1C2C3C4C5"

# =============================
# 1. Load files
# =============================

sf2 <- read_tsv(file.path(IN_DIR_SF2, "C1_BlfN_scaffold_10_wglsfs.sf2.out"),
  col_types = cols()
)

picmin <- read_tsv(file.path(IN_DIR_Pic, "pops_C1C2C3C4C5_scaffold_10.tsv"),
  col_types = cols()
)

# Keep only windows with q < 0.05
sig <- picmin %>% filter(pooled_q < 0.05)

# =============================
# 2. Annotate SF2 positions
# =============================

# Function: return q-value of PicMin window containing the SF2 location
find_q <- function(position, df) {
  hit <- df %>% filter(start <= position, end > position)
  if (nrow(hit) == 0) return(NA_real_)
  return(hit$pooled_q[1])
}

sf2$picmin_q <- sapply(sf2$location, find_q, df = sig)

sf2 <- sf2 %>%
  mutate(selection_label = ifelse(is.na(picmin_q), "neutral", "selected"))

# =============================
# 3. Save annotated table
# =============================

write_csv(sf2, "sf2_annotated_with_picmin.csv")

# =============================
# 4. Plot SweepFinder2 LR values
# =============================

cols <- c("neutral" = "grey70", "selected" = "red")

ggplot(sf2, aes(x = location, y = LR, color = selection_label)) +
  geom_point(size = 0.8, alpha = 0.9) +
  scale_color_manual(values = cols) +
  theme_minimal() +
  labs(
    title = "SweepFinder2 LR values coloured by PicMin significant windows",
    x = "Position on scaffold 10",
    y = "LR (SweepFinder2)",
    color = "Selection"
  )
