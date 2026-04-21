# ── 0. Install/load packages ─────────────────────────────────────────
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggpubr",  quietly = TRUE)) install.packages("ggpubr")

library(readxl)
library(ggplot2)
library(ggpubr)
rm(list=ls())

# ── 1. Load data ─────────────────────────────────────────────────────
df_abio <- read.csv("data/pi_abiotic.csv")
df_bio <- read.csv("data/pi_biotic.csv")

# Specify your pi column and environmental variable columns
pi_col  <- "pi"
abio_env_cols <- c("latitude", "longitude", "area", "agri_per_200", "urban_per_200", "trans", "con", "pH", "N", "P", "O", "temp")
bio_env_cols <- c("trees", "heath", "rank", "unimproved", "semi_improved", "improved", "arable", "chl_a", "phyco", "submerge", "float", "emerge")

# ── 2. Spearman correlations ──────────────────────────────────────────
results_abio <- lapply(abio_env_cols, function(var) {
  test <- cor.test(df_abio[[pi_col]], df_abio[[var]],
                   method = "spearman",
                   exact  = FALSE)   # exact=FALSE avoids ties warning
  data.frame(
    variable = var,
    rho      = round(test$estimate, 3),
    p_value  = test$p.value
  )
})
results_abio_df <- do.call(rbind, results_abio)

results_bio <- lapply(bio_env_cols, function(var) {
  test <- cor.test(df_bio[[pi_col]], df_bio[[var]],
                   method = "spearman",
                   exact  = FALSE)   # exact=FALSE avoids ties warning
  data.frame(
    variable = var,
    rho      = round(test$estimate, 3),
    p_value  = test$p.value
  )
})
results_bio_df <- do.call(rbind, results_bio)

# ── 3. Benjamini-Hochberg FDR correction ─────────────────────────────
results_abio_df$p_adj <- p.adjust(results_abio_df$p_value, method = "BH")
results_abio_df$significant <- results_abio_df$p_adj < 0.05
results_abio_df$p_value <- signif(results_abio_df$p_value, 3)
results_abio_df$p_adj   <- signif(results_abio_df$p_adj,   3)

print(results_abio_df, row.names = FALSE)
write.csv(results_abio_df, "pi_abiotic_correlations.csv", row.names = FALSE)

results_bio_df$p_adj <- p.adjust(results_bio_df$p_value, method = "BH")
results_bio_df$significant <- results_bio_df$p_adj < 0.05
results_bio_df$p_value <- signif(results_bio_df$p_value, 3)
results_bio_df$p_adj   <- signif(results_bio_df$p_adj,   3)

print(results_bio_df, row.names = FALSE)
write.csv(results_bio_df, "pi_biotic_correlations.csv", row.names = FALSE)

# ── 4. Scatter plots for significant variables ────────────────────────
sig_vars <- results_abio_df$variable[results_abio_df$significant]

# Spearman correlation: pi ~ latitude
lat_test <- cor.test(df_abio$pi, df_abio$latitude,
                     method = "spearman",
                     exact  = FALSE)

cat("rho =", round(lat_test$estimate, 3),
    " | p =", signif(lat_test$p.value, 3), "\n")

# Plot
library(ggplot2)

p <- ggplot(df_abio, aes(x = latitude, y = pi)) +
  geom_point(size = 2.5, alpha = 0.8, colour = "#1D9E75") +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#0F6E56", fill = "#9FE1CB", alpha = 0.2) +
  annotate("text",
           x = max(df_abio$latitude) * 0.999,
           y = max(df_abio$pi),
           label = paste0("rho = ", round(lat_test$estimate, 3),
                          "\np = ", signif(lat_test$p.value, 3)),
           hjust = 1, vjust = 1, size = 3.5, colour = "grey40") +
  labs(x = "Latitude", y = "Nucleotide diversity (pi)") +
  theme_bw(base_size = 12)

ggsave("plots/pi_vs_latitude.png", p, dpi=300)


# EXTRA: TajD
tajd_col  <- "tajD"
abio_env_cols <- c("latitude", "longitude", "area", "agri_per_200", "urban_per_200", "trans", "con", "pH", "N", "P", "O", "temp")
bio_env_cols <- c("trees", "heath", "rank", "unimproved", "semi_improved", "improved", "arable", "chl_a", "phyco", "submerge", "float", "emerge")

results_abio <- lapply(abio_env_cols, function(var) {
  test <- cor.test(df_abio[[tajd_col]], df_abio[[var]],
                   method = "spearman",
                   exact  = FALSE)   # exact=FALSE avoids ties warning
  data.frame(
    variable = var,
    rho      = round(test$estimate, 3),
    p_value  = test$p.value
  )
})
results_abio_df <- do.call(rbind, results_abio)

results_abio_df$p_adj <- p.adjust(results_abio_df$p_value, method = "BH")
results_abio_df$significant <- results_abio_df$p_adj < 0.05
results_abio_df$p_value <- signif(results_abio_df$p_value, 3)
results_abio_df$p_adj   <- signif(results_abio_df$p_adj,   3)

print(results_abio_df, row.names = FALSE)
write.csv(results_abio_df, "tajD_abiotic_correlations.csv", row.names = FALSE)
