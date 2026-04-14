### ============================================
### Load libraries
### ============================================
library(tidyverse)
library(data.table)
library(dplyr)


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
tajD_df <- read_pixy_stat("tajima_d", "pixy_stats")
theta_df <- read_pixy_stat("watterson_theta", "pixy_stats")

# Remove windows with NA or no SNPs
pi_df <- pi_df %>% filter(!is.na(avg_pi), no_sites > 0)
tajD_df <- tajD_df %>% filter(!is.na(tajima_d), no_sites > 0)
theta_df <- theta_df %>% filter(!is.na(raw_watterson_theta), no_sites > 0)

summary(pi_df)
summary(tajD_df)
summary(theta_df)


### ============================================
### Genome‑wide averages per population & merge
### ============================================
genomewide_summary <- list(
  
  # pi
  pi_df %>%
    group_by(pop) %>%
    summarise(
      pi_mean_unweighted = mean(avg_pi, na.rm = TRUE),
      pi_mean_weighted   = sum(avg_pi * no_sites, na.rm = TRUE) /
        sum(no_sites, na.rm = TRUE),
      n_windows_pi       = n(),
      .groups = "drop"),
  # Watterson's theta
  theta_df %>%
    group_by(pop) %>%
    summarise(
      thetaW_mean_unweighted = mean(avg_watterson_theta, na.rm = TRUE),
      thetaW_mean_weighted   = sum(avg_watterson_theta * no_sites, na.rm = TRUE) /
        sum(no_sites, na.rm = TRUE),
      n_windows_thetaW       = n(),
      .groups = "drop"),
  # Tajima's D
  tajD_df %>%
    group_by(pop) %>%
    summarise(
      tajD_mean_unweighted = mean(tajima_d, na.rm = TRUE),
      tajD_mean_weighted   = sum(tajima_d * no_sites, na.rm = TRUE) /
        sum(no_sites, na.rm = TRUE),
      n_windows_tajD       = n(),
      .groups = "drop")
) %>%
  reduce(left_join, by = "pop")


### ============================================
### Plot weighted pi and theta, unweighted tajima's d
### ============================================

# all 3 metrics
plot_df <- genomewide_summary %>%
  select(
    pop,
    pi_mean_weighted,
    thetaW_mean_weighted,
    tajD_mean_unweighted
  ) %>%
  rename(
    pi = pi_mean_weighted,
    thetaW = thetaW_mean_weighted,
    tajimaD = tajD_mean_unweighted
  ) %>%
  pivot_longer(
    cols = c(pi, thetaW, tajimaD),
    names_to = "metric",
    values_to = "value"
  )

ggplot(plot_df, aes(x = pop, y = value)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ metric, scales = "free_y",
             labeller = as_labeller(c(
               pi = "Weighted π",
               thetaW = "Weighted θW",
               tajimaD = "Unweighted Tajima’s D"
             ))) +
  labs(
    x = "Population",
    y = "Genome-wide value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

### only tajima's D
tajd_plot_df <- genomewide_summary %>%
  select(pop, tajimaD = tajD_mean_unweighted)

# ordered barplot
p_tajD<-ggplot(
  tajd_plot_df,
  aes(x = reorder(pop, tajimaD), y = tajimaD)
) +
  geom_col(aes(fill = tajimaD > 0), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("firebrick", "steelblue")) +
  labs(
    x = "Population",
    y = "Mean Tajima’s D (unweighted)"
  ) +
  theme_bw()

ggsave("tajimasD_ordered_barplot.png", plot = p_tajD, dpi = 300)
### ============================================
### Check for correlations
### ============================================

corr_df <- genomewide_summary %>%
  select(
    pop,
    pi = pi_mean_weighted,
    thetaW = thetaW_mean_weighted,
    tajimaD = tajD_mean_unweighted
  )
# π vs θ: expect high correlation
# Pearson → linear relationship
# Spearman → rank‑based (robust to outliers)
cor.test(corr_df$pi, corr_df$thetaW, method = "pearson")
cor.test(corr_df$pi, corr_df$thetaW, method = "spearman")

# π vs tajD: Spearman because Tajima’s D is not linear and not normally distributed
cor.test(corr_df$pi, corr_df$tajimaD, method = "spearman")

# θ vs tajD: Spearman
cor.test(corr_df$thetaW, corr_df$tajimaD, method = "spearman")

# correlation matrix
cor_mtx <-cor(corr_df %>% select(-pop), method = "spearman")


### ============================================
### Scatterplot π vs θ
### ============================================
plot_df <- genomewide_summary %>%
  select(
    pop,
    pi = pi_mean_weighted,
    thetaW = thetaW_mean_weighted
  )
# plot a 
ggplot(plot_df, aes(x = thetaW, y = pi, label = pop)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(nudge_y = 0.00005, size = 3) +
  labs(
    x = "Watterson’s θ (weighted)",
    y = "Nucleotide diversity π (weighted)"
  ) +
  theme_bw()

# plot b -better
p_pitheta <- ggplot(plot_df, aes(x = thetaW, y = pi, label = pop)) +
  geom_point(size = 3, color = "steelblue") +
  
  ## θ = π reference line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  ## population labels (close to points)
  geom_text(
    vjust = -1,   # small vertical offset
    hjust = 0.5,
    size = 3
  ) +
  
  ## optional linear regression
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  
  labs(
    x = "Watterson’s θ (weighted)",
    y = "Nucleotide diversity π (weighted)"
  ) +
  
  theme_bw()

ggsave("pi_vs_theta_scatterplot.png", plot = p_pitheta, dpi = 300)

