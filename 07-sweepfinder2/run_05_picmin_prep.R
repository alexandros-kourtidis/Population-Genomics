rm(list = ls())
#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)

setwd("C:/Users/u0167113/Documents/PhD research/popgen analyses/07. SweepFinder2/polarised/")

# --------------------------
# USER SETTINGS
# --------------------------

populations <- c(
  "A1_AnOudin", "A2_BuSN", "A3_Mech", "A4_PL15_YEL", "A5_GenM",
  "B1_BKN1", "B2_OM2", "B3_ZW", "B4_OHZ", "B5_DA2",
  "C1_BlfN", "C2_MO", "C3_Ter1", "C4_BW_48630", "C5_BW_36962",
  "D1_CBOO6", "D2_LRV", "D3_BKLE5", "D4_BW_62256", "D5_BW_22050"
)

base_dir <- "sf2_out"
window_size <- 10000   # 10 kb


# ================================================
# STEP 1: COMBINE SCAFFOLD-LEVEL SF2 OUTPUTS
# ================================================

for (pop in populations) {
  
  message("Combining SF2 outputs for: ", pop)
  
  pop_dir <- file.path(base_dir, pop)
  
  files <- list.files(pop_dir,
                      pattern = paste0(pop, ".*sf2.out$"),
                      full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No SF2 files found for ", pop)
    next
  }
  
  # Read and combine files
  dat_list <- lapply(files, function(f) {
    scaf <- str_extract(basename(f), "scaffold_[0-9]+")
    df <- read_table2(f, col_names = TRUE, comment = "#")
    names(df) <- tolower(names(df))  # normalize column names
    df$scaffold <- scaf
    df
  })
  
  combined <- bind_rows(dat_list)
  
  out_file <- file.path(pop_dir, paste0("combined_", pop, ".txt"))
  write.table(combined, out_file, row.names = FALSE, quote = FALSE, sep = "\t")
  
  message("Written: ", out_file)
}

# ========================================
# STEP 2: WINDOW AGGREGATION
# ========================================

for (pop in populations) {
  
  message("Windowing for: ", pop)
  
  pop_dir <- file.path(base_dir, pop)
  combined_file <- file.path(pop_dir, paste0("combined_", pop, ".txt"))
  
  if (!file.exists(combined_file)) {
    warning("Missing combined file for ", pop, " — skipping.")
    next
  }
  
  combined <- read.table(combined_file, header = TRUE)
  
  # Expect columns: location, lr
  names(combined) <- tolower(names(combined))
  
  if (!("location" %in% names(combined)) |
      !("lr" %in% names(combined))) {
    stop("combined file for ", pop, " does not contain location or lr columns!")
  }
  
  # Add window index
  combined <- combined %>%
    mutate(window = floor(location / window_size))
  
  # Compute mean LR per window
  window_summary <- combined %>%
    group_by(scaffold, window) %>%
    summarise(mean_LR = mean(lr, na.rm = TRUE),
              n_sites = n(),
              .groups = "drop") %>%
    arrange(scaffold, window)
  
  out_file <- file.path(pop_dir, paste0(pop, "_10kb_windows_meanLR.txt"))
  write.table(window_summary, out_file,
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  message("Written: ", out_file)
}

# ========================================
# STEP 3: READ ALL WINDOW SUMMARY FILES
# ========================================

all_pops_data <- list()

for (pop in populations) {
  
  pop_dir <- file.path(base_dir, pop)
  infile  <- file.path(pop_dir, paste0(pop, "_10kb_windows_meanLR.txt"))
  
  if (!file.exists(infile)) {
    warning("Missing file: ", infile)
    next
  }
  
  df <- read.table(infile, header = TRUE, sep = "\t") %>%
    mutate(pop = pop)
  
  all_pops_data[[pop]] <- df
}

# Combine all into one long data frame
long_df <- bind_rows(all_pops_data)

# ================================
# STEP 4: CREATE PICMIN-READY INPUT
# ================================
# We want: 

library(data.table)
library(dplyr)
library(tidyr)

WIN_SIZE <- window_size 

picmin_input <- long_df %>%
  select(scaffold, window, pop, mean_LR) %>%
  pivot_wider(
    id_cols  = c(scaffold, window),
    names_from  = pop,
    values_from = mean_LR
  ) %>%
  arrange(scaffold, window) %>%
  mutate(
    start     = window * WIN_SIZE,
    end       = window * WIN_SIZE + WIN_SIZE,
    window_id = paste0(scaffold, "_", start)
  )

# Reorder columns: window_id first (preferred for PicMin workflow)
picmin_input <- picmin_input %>%
  relocate(window_id, .before = scaffold)

# Create output directory
out_dir <- "picmin"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- "PicMin_input_meanLR_10kb_windows_woContam.txt"
outfile <- file.path(out_dir, out_file)

write.table(
  picmin_input,
  outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
