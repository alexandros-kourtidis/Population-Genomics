rm(list=ls())
# ================================
# Load libraries
# ================================
library(tidyverse)

# ================================
# Set directory containing .het files
# ================================
het_dir <- "het_stats"   # change if needed

# List all het files
het_files <- list.files(, pattern = "\\.het$", full.names = TRUE)

# Output file
outfile <- "F_values_all_populations.png"

# ================================
# Read and combine all het files
# ================================
het_data <- het_files %>%
  map_df(function(file) {
    
    # Extract population from filename
    # Example filename: A2_BuSN_CN.filtered.gtonly.het_per_sample.het
    pop <- basename(file) %>% str_replace(".filtered.*", "")
    
    df <- read.table(file, header = TRUE)
    
    df$Population <- pop
    
    # Clean sample name (remove ".bam" and directories)
    df$Sample <- df$INDV %>%
      basename() %>% 
      str_replace(".bam$", "")
    
    df
  })

# ================================
# Barplot of F values for each population
# Creates one PNG per population
# ================================
outdir <- "plots_F_per_population"
dir.create(outdir, showWarnings = FALSE)

for (pop in unique(het_data$Population)) {
  
  sub <- het_data %>% filter(Population == pop)
  
  p <- ggplot(sub, aes(x = Sample, y = F)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    ggtitle(paste("Inbreeding coefficient (F) per sample -", pop)) +
    ylab("F value") +
    xlab("Sample")
  
  ggsave(filename = file.path(outdir, paste0(pop, "_F_barplot.png")),
         plot = p, width = 10, height = 6, dpi = 300)
}

cat("Finished! Plots saved in:", outdir, "\n")


# ======================================
# Combined faceted barplot with red asterisks
# ======================================
# Define the samples to annotate
special_samples <- c(
  "A1_AnOudin_C14_aligned.sorted.nd.bam",
  "D2_LRV_C32_aligned.sorted.nd.bam"
)

# Clean them the same way as in the dataframe
special_samples_clean <- special_samples %>%
  basename() %>%
  str_replace(".bam$", "")

# Compute y-position for asterisks (slightly below bars)
# You can also put them ABOVE the bars if you prefer.
asterisk_df <- het_data %>%
  filter(Sample %in% special_samples_clean) %>%
  mutate(
    y_pos = pmin(F, 0) - 0.07   # place slightly BELOW bar
    # If you want ABOVE, use: y_pos = pmax(F, 0) + 0.02
  )

# --------------Plot---------------------
# Order populations so layout stays consistent
het_data$Population <- factor(het_data$Population,
                              levels = unique(het_data$Population))

# Number of rows you want:
n_rows <- 3   # change to 2 if you want 2 rows instead

p <- ggplot(het_data, aes(x = Sample, y = F)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(
    data = asterisk_df,
    aes(x = Sample, y = y_pos, label = "*"),
    color = "red",
    size = 8,
    fontface = "bold"
  ) +
  facet_wrap(~ Population, scales = "free_x", nrow = n_rows) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  ylab("Inbreeding coefficient (F)") +
  xlab("Sample") +
  ggtitle("F values (heterozygosity / inbreeding) across all populations")

# Save high-resolution figure
ggsave(outfile, p, width = 20, height = 10, dpi = 300)

cat("Plot saved as:", outfile, "\n")
