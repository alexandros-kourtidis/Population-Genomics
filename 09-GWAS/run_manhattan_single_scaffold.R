
# install.packages("viridis")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

# --- Parameters (edit as needed) ---
gwas_file      <- "batch1234_plink_gemma_lmm_results_nocov.assoc.txt"  # GEMMA output
scaffold_file  <- "scaffold_lengths.txt"                                # Scaffold lengths (two columns: name, length_bp)
output_png     <- "manhattan_plot_scaffold6.png"                        # Output image
downsample     <- 1e5                                                  # e.g., 1e5, or NULL to plot all
scaffold_to_plot <- 6                                                   # can be 6, "6", or "scaffold_6"

# --- Helper: normalize scaffold name to match data ---
normalize_scaffold <- function(x, scaffold_names) {
  # If user gave a number (6), accept "6" and "scaffold_6"
  if (is.numeric(x)) {
    candidates <- c(as.character(x), paste0("scaffold_", x))
  } else {
    # If user gave "6" or "scaffold_6", consider both forms
    sx <- as.character(x)
    num <- sub("^scaffold_", "", sx)
    candidates <- unique(c(sx, num, paste0("scaffold_", num)))
  }
  # Pick the first candidate that exists in the dataset
  match <- candidates[candidates %in% scaffold_names]
  if (length(match) == 0) {
    stop("Requested scaffold not found. Tried: ", paste(candidates, collapse = ", "))
  }
  match[1]
}

# --- Load scaffold lengths ---
scaffolds <- fread(scaffold_file)[, 1:2]
setnames(scaffolds, c("scaffold", "length"))

# Make sure length is numeric
scaffolds <- scaffolds %>%
  mutate(length = as.numeric(length))

# Normalize requested scaffold name to what's present
scaffold_name <- normalize_scaffold(scaffold_to_plot, scaffolds$scaffold)

# Filter to the single scaffold and compute local positions (start at 0)
scaffolds_single <- scaffolds %>%
  filter(scaffold == scaffold_name) %>%
  mutate(chromStart = 0,
         chromMid   = length / 2)

# --- Load GWAS results ---
gwas <- fread(gwas_file, select = c("chr", "ps", "p_wald"))
gwas[, ps := as.numeric(ps)]
gwas[, p_wald := as.numeric(p_wald)]

# Filter GWAS to the chosen scaffold and compute positions
gwas_single <- gwas %>%
  filter(chr == scaffold_name) %>%
  left_join(scaffolds_single, by = c("chr" = "scaffold")) %>%
  mutate(pos_cum = ps + chromStart)

# Safety check
if (nrow(gwas_single) == 0) {
  stop("No SNPs found on scaffold: ", scaffold_name,
       ". Check that 'chr' in GWAS matches scaffold names in ", scaffold_file, ".")
}

# --- Significance thresholds using the filtered set ---
num_tests <- nrow(gwas_single)
bonferroni_threshold <- 0.05 / num_tests
bonferroni_log <- -log10(bonferroni_threshold)
nominal_log <- -log10(0.05)

# --- Optional downsampling (after filtering) ---
if (!is.null(downsample) && nrow(gwas_single) > downsample) {
  set.seed(123)
  gwas_plot <- gwas_single[sample(.N, downsample)][order(pos_cum)]
} else {
  gwas_plot <- gwas_single[order(pos_cum)]
}

# --- Plot to PNG (single scaffold) ---
png(output_png, width = 2000, height = 1200, res = 300)

ggplot(gwas_plot, aes(x = pos_cum, y = -log10(p_wald))) +
  geom_point(color = "darkgray", alpha = 0.6, size = 0.6) +
  scale_x_continuous(
    breaks = scaffolds_single$chromMid,
    labels = scaffolds_single$scaffold
  ) +
  geom_hline(yintercept = bonferroni_log, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = nominal_log, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste("Manhattan Plot (single scaffold):", scaffold_name),
    x = "Scaffold",
    y = expression(-log[10](p))
  ) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title    = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 14)
  )

dev.off()

