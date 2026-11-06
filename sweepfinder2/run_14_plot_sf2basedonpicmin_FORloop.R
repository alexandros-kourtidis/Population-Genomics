rm(list = ls())
library(dplyr)
library(readr)
library(ggplot2)
library(fuzzyjoin)
library(purrr)
library(stringr)

setwd("C:/Users/u0167113/Documents/PhD research/popgen analyses/07. SweepFinder2/polarised")

# ===================================================
# 0. Variables and filters
# ===================================================

SCAFFOLD        <- "scaffold_10"
Q_THRESHOLD     <- 0.05
N_EST_THRESHOLD <- 2       # <-- set to NA if you don’t want an n_est filter
SUFFIX          <- paste0("q", Q_THRESHOLD, if (!is.na(N_EST_THRESHOLD)) paste0("_n", N_EST_THRESHOLD) else "")

IN_DIR_SF2 <- "sf2_out"
IN_DIR_PIC <- "picmin/results/pops_C1C2C3C4C5"

pops <- c(
  #"A1_AnOudin", "A2_BuSN", "A3_Mech", "A4_PL15_YEL", "A5_GenM",
  #"B1_BKN1", "B2_OM2", "B3_ZW", "B4_OHZ", "B5_DA2",
  "C1_BlfN", "C2_MO", "C3_Ter1", "C4_BW_48630", "C5_BW_36962"
  #,
  #"D1_CBOO6", "D2_LRV", "D3_BKLE5", "D4_BW_62256", "D5_BW_22050"
)

# ===================================================
# 1. Load PicMin file once
# ===================================================

picmin <- read_tsv(
  file.path(IN_DIR_PIC, paste0("pops_C1C2C3C4C5_", SCAFFOLD, ".tsv")),
  col_types = cols()
) %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end),
    pooled_q = as.numeric(pooled_q),
    n_est    = as.numeric(n_est)
  )

# ===================================================
# 2) Build a MASTER annotated table
#    - join SF2 positions to ALL PicMin windows
#    - if multiple windows overlap a position, keep:
#        picmin_q_min = min(pooled_q)
#        n_est_max    = max(n_est)
# ===================================================

# Function to
build_master_for_pop <- function(pop) {
  sf2_file <- file.path(IN_DIR_SF2, pop, paste0(pop, "_", SCAFFOLD, "_wglsfs.sf2.out"))
  message("Joining PicMin → ", pop)
  sf2 <- read_tsv(sf2_file, col_types = cols())
  
  # range join: position ∈ [start, end]
  hits <- fuzzy_inner_join(
    sf2, picmin,
    by = c("location" = "start", "location" = "end"),
    match_fun = list(`>=`, `<=`)
  )
  
  # If a location overlaps multiple PicMin windows, summarise to a single row per location
  picmin_summarised <- hits %>%
    group_by(location) %>%
    summarise(
      picmin_q_min = suppressWarnings(min(pooled_q, na.rm = TRUE)),
      n_est_max    = suppressWarnings(max(n_est, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      picmin_q_min = ifelse(is.infinite(picmin_q_min), NA_real_, picmin_q_min),
      n_est_max    = ifelse(is.infinite(n_est_max),    NA_real_, n_est_max)
    )
  
  # Left join summaries back to SF2; keep the original LR & other SF2 columns
  annotated <- sf2 %>%
    left_join(picmin_summarised, by = "location") %>%
    mutate(population = pop)
  
  annotated
}

# Build master once for all populations
master <- map_dfr(pops, build_master_for_pop)

# Optionally save an unfiltered master per population and/or a combined master
# (This file can be reused for any thresholds without rejoining)
write_csv(master, paste0("sf2_master_", SCAFFOLD, ".csv"))


# =========================================================
# 3) Thresholds applied only for classification/plot/export
# =========================================================
apply_thresholds_and_label <- function(df, q_thr, n_thr) {
  df %>%
    mutate(
      passes_q   = ifelse(!is.na(q_thr), !is.na(picmin_q_min) & picmin_q_min <= q_thr, !is.na(picmin_q_min)),
      passes_n   = ifelse(!is.na(n_thr), !is.na(n_est_max)    & n_est_max    >= n_thr, TRUE),
      selected   = passes_q & passes_n,
      selection_label = ifelse(selected, "selected", "neutral")
    )
}

annotated <- apply_thresholds_and_label(master, Q_THRESHOLD, N_EST_THRESHOLD)

# Save per-population annotated tables WITH labels (but built from the master)
invisible(
  annotated %>%
    group_split(population) %>%
    walk(function(dfp) {
      pop <- unique(dfp$population)
      out <- paste0("sf2_annotated_", SCAFFOLD, "_", pop, "_", SUFFIX, ".csv")
      write_csv(dfp, out)
    })
)

# =========================================================
# 4) Plot
# =========================================================
cols <- c("neutral" = "grey70", "selected" = "red")

p <- ggplot(annotated, aes(x = location, y = LR, color = selection_label)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(values = cols) +
  facet_wrap(~ population, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  labs(
    title = paste0(
      "SweepFinder2 LR coloured by PicMin windows (", SCAFFOLD, ")\n",
      "Filters: q ≤ ", Q_THRESHOLD,
      if (!is.na(N_EST_THRESHOLD)) paste0(", n_est ≥ ", N_EST_THRESHOLD) else " (no n_est filter)"
    ),
    x = paste("Position on", SCAFFOLD),
    y = "LR"
  )

p

ggsave(paste0("SF2_PicMin_plot_", SCAFFOLD, "_", SUFFIX, ".png"), p, width = 10, height = 14, dpi = 300)


