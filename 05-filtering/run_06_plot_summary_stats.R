#!/usr/bin/env Rscript
# ============================================================
# compare_vcf_stats_combined.R
#
# Compares bcftools stats output (SN / TSTV / PSC blocks) across
# filtering steps for either invariants-included
# run, or SNPs-only run. Run this script twice, changing only
# `inv_status` (and `out_dir`, which is derived from it) between runs.
#
# Everything you might need to change lives in the CONFIG block
# below -- directories, filenames, filtering steps, and output
# plot/table names. Nothing else in the script should need editing
# unless your bcftools stats file layout itself changes.
# ============================================================
rm(list = ls())

library(tidyverse)

## ============================================================
## ---- CONFIG -- edit everything here ------------------------
## ============================================================

stats_dir <- "summaries"
prefix    <- "batch1234_merged1_10"

# Define every stage you want to compare, as (inv_status, stage) pairs.
stage_order <- c(
  "raw",
  "inv.DPMQ_DP_MISS_noMitoMiss2",
  "snp.QUALDPMQ_DPGQ_MISS_noMitoMiss2",
  "snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss2",
  "invsnp.filtered_sorted2"
)

out_dir <- file.path("plots", "06_combined_aftercleannames")

csv_sn_summary        <- "SN_summary_by_step.csv"
csv_tstv_summary      <- "TsTv_by_step.csv"
csv_persample_missing <- "per_sample_missingness_by_step.csv"
csv_overall_missing   <- "overall_missingness_by_step.csv"
csv_imiss_summary     <- "imiss_summary_by_stage.csv"
csv_worst_samples     <- "imiss_worst_samples.csv"

png_variant_counts <- "01_variant_counts.png"
png_tstv_ratio     <- "02_tstv_ratio.png"
png_depth_boxplot  <- "03_depth_boxplot.png"
png_snp_counts     <- "04_snp_counts.png"
png_imiss_points   <- "05_imiss_per_sample_points.png"
png_imiss_slope    <- "06_imiss_slope_by_stage.png"
png_imiss_heatmap  <- "07_imiss_heatmap.png"
png_imiss_violin   <- "08_imiss_violin_by_stage.png"

keep_keys <- c("number of records", "number of SNPs", "number of indels",
               "number of multiallelic sites", "number of multiallelic SNP sites",
               "number of no-ALTs", "number of MNPs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## ---- File list (built from stage_order) -----------------------
## ============================================================
# Expected filename: <prefix>.<stage>.summary.txt

stats_files <- file.path(stats_dir, paste0(prefix, ".", stage_order, ".summary.txt"))
names(stats_files) <- stage_order

missing <- stats_files[!file.exists(stats_files)]
if (length(missing) > 0) {
  stop("These files were not found:\n", paste(missing, collapse = "\n"))
}

## ============================================================
## ---- Parsers for the relevant bcftools stats blocks ----------
## ============================================================

parse_sn <- function(lines) {
  sn_lines <- lines[grepl("^SN\t", lines)]
  if (length(sn_lines) == 0) return(tibble())
  do.call(rbind, lapply(sn_lines, function(l) {
    f <- str_split(l, "\t")[[1]]
    tibble(key = str_remove(f[3], ":$"), value = as.numeric(f[4]))
  }))
}

parse_tstv <- function(lines) {
  l <- lines[grepl("^TSTV\t", lines)]
  if (length(l) == 0) return(tibble(ts = NA_real_, tv = NA_real_))
  f <- str_split(l[1], "\t")[[1]]
  tibble(ts = as.numeric(f[3]), tv = as.numeric(f[4]))
}

parse_psc <- function(lines) {
  l <- lines[grepl("^PSC\t", lines)]
  if (length(l) == 0) return(tibble())
  do.call(rbind, lapply(l, function(x) {
    f <- str_split(x, "\t")[[1]]
    tibble(
      sample         = f[3],
      nRefHom        = as.numeric(f[4]),
      nNonRefHom     = as.numeric(f[5]),
      nHets          = as.numeric(f[6]),
      nTransitions   = as.numeric(f[7]),
      nTransversions = as.numeric(f[8]),
      nIndels        = as.numeric(f[9]),
      avgDepth       = as.numeric(f[10]),
      nSingletons    = as.numeric(f[11]),
      nMissing       = as.numeric(f[14])
    )
  }))
}

## ============================================================
## ---- Load all stats files (one per stage) --------------
## ============================================================

sn_all   <- list()
tstv_all <- list()
psc_all  <- list()

for (stg in stage_order) {
  f     <- stats_files[[stg]]
  lines <- readLines(f, warn = FALSE)
  
  sn <- parse_sn(lines)
  if (nrow(sn) > 0) {
    sn$stage <- stg
    sn_all[[stg]] <- sn
  }
  
  tv <- parse_tstv(lines)
  tv$stage <- stg
  tstv_all[[stg]] <- tv
  
  psc <- parse_psc(lines)
  if (nrow(psc) > 0) {
    psc$stage <- stg
    psc_all[[stg]] <- psc
  }
}

sn_df   <- bind_rows(sn_all)   %>% mutate(stage = factor(stage, levels = stage_order))
tstv_df <- bind_rows(tstv_all) %>% mutate(stage = factor(stage, levels = stage_order))
psc_df  <- bind_rows(psc_all)  %>% mutate(stage = factor(stage, levels = stage_order))

## ============================================================
## ---- Reshape (one file per stage -- nothing to sum) ----
## ============================================================

sn_agg <- sn_df %>% distinct(stage, key, value)
sn_agg_plot <- sn_agg %>% filter(key %in% keep_keys)

tstv_agg <- tstv_df %>%
  distinct(stage, ts, tv) %>%
  mutate(ts_tv = ts / tv)

psc_agg_sample <- psc_df %>%
  mutate(
    total_calls   = nRefHom + nNonRefHom + nHets + nMissing,
    missing_frac  = nMissing / total_calls,
    het_hom_ratio = nHets / pmax(nNonRefHom, 1)
  ) %>%
  select(stage, sample, nMissing, total_calls, avgDepth,
         nHets, nNonRefHom, missing_frac, het_hom_ratio)

psc_agg_overall <- psc_agg_sample %>%
  group_by(stage) %>%
  summarise(
    mean_missing_frac = mean(missing_frac, na.rm = TRUE),
    sd_missing_frac   = sd(missing_frac, na.rm = TRUE),
    mean_depth        = mean(avgDepth, na.rm = TRUE),
    .groups = "drop"
  )

# imiss-equivalent table: one row per sample per filter step, with a
# fractional-missingness column -- plays the role of a real .imiss file.
all_data <- psc_agg_sample %>%
  transmute(INDV = sample, stage, F_MISS = missing_frac)

## ============================================================
## ---- Write tables ---------------------------------------------
## ============================================================

write_csv(sn_agg,          file.path(out_dir, csv_sn_summary))
write_csv(tstv_agg,        file.path(out_dir, csv_tstv_summary))
write_csv(psc_agg_sample,  file.path(out_dir, csv_persample_missing))
write_csv(psc_agg_overall, file.path(out_dir, csv_overall_missing))

## ============================================================
## ---- bcftools-derived plots -------------------------------------
## ============================================================

p_counts <- ggplot(sn_agg_plot, aes(x = stage, y = value, fill = stage)) +
  geom_col() +
  geom_text(aes(label = scales::comma(value)), vjust = 1.3, color = "white", size = 2.3) +
  facet_wrap(~ key, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Variant counts across filtering steps (log scale)",
       x = "Filtering step", y = "Count (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, png_variant_counts), p_counts, width = 12, height = 6, dpi = 200)

p_tstv <- ggplot(tstv_agg, aes(x = stage, y = ts_tv, group = 1)) +
  geom_line(alpha = 0.4, color = "grey50") +
  geom_point(aes(color = stage), size = 3) +
  labs(title = "Ts/Tv ratio across filtering steps",
       x = "Filtering step", y = "Ts/Tv") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, png_tstv_ratio), p_tstv, width = 7, height = 5, dpi = 150)

p_depth <- ggplot(psc_agg_sample, aes(x = stage, y = avgDepth, fill = stage)) +
  geom_boxplot() +
  labs(title = "Per-sample average depth across filtering steps",
       x = "Filtering step", y = "Average depth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, png_depth_boxplot), p_depth, width = 7, height = 5, dpi = 150)

snp_counts <- sn_agg %>%
  filter(key == "number of SNPs") %>%
  mutate(stage = factor(stage, levels = rev(stage_order)))

p_snp_counts <- ggplot(snp_counts, aes(y = stage, x = value, fill = stage)) +
  geom_col() +
  geom_text(aes(label = scales::comma(value)), hjust = 1.1, color = "white", size = 3.2) +
  scale_fill_manual(values = setNames(scales::hue_pal()(length(stage_order)), rev(stage_order))) +
  scale_x_continuous(labels = scales::comma, position = "top", expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Number of SNPs across filtering steps",
       x = "Number of SNPs", y = "Filtering step") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, png_snp_counts), p_snp_counts, width = 7, height = 5, dpi = 150)

## ============================================================
## ---- imiss-equivalent plots ------------------------------------
## ============================================================

n_samp <- n_distinct(all_data$INDV)

p_points <- ggplot(all_data, aes(x = reorder(INDV, F_MISS), y = F_MISS, color = stage)) +
  geom_point(size = 2) +
  coord_flip() +
  labs(
    title = "Per-sample missingness across filtering steps",
    x = "Sample",
    y = "Fraction missing (F_MISS)",
    color = "Filtering step"
  ) +
  theme_minimal()
ggsave(file.path(out_dir, png_imiss_points), p_points,
       width = 8, height = max(4, 0.10 * n_samp), dpi = 300, limitsize = FALSE)

p_slope <- ggplot(all_data, aes(x = stage, y = F_MISS, group = INDV)) +
  geom_line(alpha = 0.4, color = "grey50") +
  geom_point(aes(color = stage), size = 2) +
  labs(
    title = "Missingness per sample across filtering steps",
    x = "Filtering step",
    y = "Fraction missing (F_MISS)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, png_imiss_slope), p_slope, width = 8, height = 5, dpi = 300)

p_heatmap <- ggplot(all_data, aes(x = stage, y = reorder(INDV, F_MISS), fill = F_MISS)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(
    title = "Missingness heatmap: samples x filtering steps",
    x = "Filtering step",
    y = "Sample",
    fill = "F_MISS"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, png_imiss_heatmap), p_heatmap,
       width = 8, height = max(4, 0.10 * n_samp), dpi = 300)

p_violin <- ggplot(all_data, aes(x = stage, y = F_MISS, fill = stage)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  labs(
    title = "Distribution of per-sample missingness by filtering step",
    x = "Filtering step",
    y = "Fraction missing (F_MISS)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, png_imiss_violin), p_violin, width = 8, height = 5, dpi = 300)

## ============================================================
## ---- Summary tables ---------------------------------------------
## ============================================================

imiss_summary_table <- all_data %>%
  group_by(stage) %>%
  summarise(
    n_samples          = n(),
    mean_fmiss         = mean(F_MISS),
    median_fmiss       = median(F_MISS),
    max_fmiss          = max(F_MISS),
    n_samples_gt_10pct = sum(F_MISS > 0.1),
    .groups = "drop"
  )
print(imiss_summary_table)
write_csv(imiss_summary_table, file.path(out_dir, csv_imiss_summary))

# Samples whose missingness is consistently high across ALL filtering
# steps (candidates for removal regardless of SNP filtering thresholds)
worst_samples <- all_data %>%
  group_by(INDV) %>%
  summarise(mean_fmiss_all_steps = mean(F_MISS), .groups = "drop") %>%
  arrange(desc(mean_fmiss_all_steps))
print(head(worst_samples, 20))
write_csv(worst_samples, file.path(out_dir, csv_worst_samples))

