#install.packages("viridis")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

# --- Parameters (edit as needed) ---
gwas_file <- "batch1234_plink_gemma_lmm_results_nocov.assoc.txt"  # GEMMA output
scaffold_file <- "scaffold_lengths.txt"                         # Scaffold lengths from VCF header
output_png <- "manhattan_plot.png"                           # Output image
downsample <- 1e5                                        # number of SNPs to plot; set NULL to plot all

# --- Load scaffold lengths ---
scaffolds <- fread(scaffold_file)[, 1:2]
setnames(scaffolds, c("scaffold", "length"))

scaffolds <- scaffolds %>%
  mutate(
    scaffold_num = as.numeric(gsub("scaffold_", "", scaffold)),
    length = as.numeric(length)
  ) %>%
  arrange(scaffold_num) %>%
  mutate(
    chromStart = lag(cumsum(length), default = 0),
    chromMid = chromStart + length / 2
  )

# --- Load GWAS results ---
gwas <- fread(gwas_file, select = c("chr", "ps", "p_wald"))
gwas[, ps := as.numeric(ps)]
gwas[, p_wald := as.numeric(p_wald)]

# --- Compute cumulative positions ---
gwas <- gwas %>%
  left_join(scaffolds, by = c("chr" = "scaffold")) %>%
  mutate(pos_cum = ps + chromStart)

# --- Significance thresholds ---
num_tests <- nrow(gwas)
bonferroni_threshold <- 0.05 / num_tests
bonferroni_log <- -log10(bonferroni_threshold)
nominal_log <- -log10(0.05)

# --- Optional downsampling ---
if(!is.null(downsample) & nrow(gwas) > downsample){
  set.seed(123)
  gwas_plot <- gwas[sample(.N, downsample)]
} else {
  gwas_plot <- gwas
}

# --- Colors for scaffolds ---
chrom_colors <- rep(c("grey40", "grey70"), length.out = nrow(scaffolds))

# --- Plot to PNG ---
png(output_png, width = 3000, height = 1500, res = 300)  # high-res image

ggplot(gwas_plot, aes(x = pos_cum, y = -log10(p_wald))) +
  geom_point(aes(color = factor(chr)), alpha = 0.6, size = 0.5) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaffolds$chromMid,
    labels = scaffolds$scaffold
  ) +
  geom_hline(yintercept = bonferroni_log, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = nominal_log, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Manhattan Plot: Daphnia magna GWAS",
    x = "Scaffold",
    y = expression(-log[10](p))
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  )
dev.off()


