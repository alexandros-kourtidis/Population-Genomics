# =============================================================================
#   Haplotype diversity estimators per population
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

rm(list=ls())

# ── 1. Core diversity functions ───────────────────────────────────────────────

# Haplotype richness (number of unique haplotypes present)
hap_richness <- function(counts) sum(counts > 0)

# Haplotype diversity (Nei 1987) — equivalent to expected heterozygosity
# Hd = (n / n-1) * (1 - sum(pi^2))
# where pi = frequency of haplotype i
hap_diversity <- function(counts) {
  n  <- sum(counts)
  if (n < 2) return(NA)
  p  <- counts[counts > 0] / n
  (n / (n - 1)) * (1 - sum(p^2))
}

# Simpson's index (probability that two random draws are different haplotypes)
simpson <- function(counts) {
  n <- sum(counts)
  if (n < 2) return(NA)
  1 - sum(counts * (counts - 1)) / (n * (n - 1))
}

# ── 2. Apply to each population ───────────────────────────────────────────────
# df has rows = haplotypes, columns = populations (plus Haplotype column)
df <- read.csv("pop_mitohaplotypes_matrix.csv")
populations <- colnames(df[2:21])

# Minimum sample size across populations (for rarefaction)
pop_totals <- colSums(df[, populations])
n_min      <- min(pop_totals[pop_totals > 0])
message("Rarefaction target sample size (n_min): ", n_min)

diversity_df <- do.call(rbind, lapply(populations, function(pop) {
  counts <- df[[pop]]
  data.frame(
    Population       = pop,
    N                = sum(counts),
    Richness         = hap_richness(counts),
    Hap_Diversity    = round(hap_diversity(counts), 4),
    Simpson          = round(simpson(counts),         4)
  )
}))

rownames(diversity_df) <- NULL
print(diversity_df)
write.csv(diversity_df, "mitohaplotype_diversity.csv", row.names = FALSE)

# ── 3. Plot all estimators as a bar chart ─────────────────────────────────────
diversity_long <- diversity_df |>
  pivot_longer(
    cols      = c(Richness, Hap_Diversity, Simpson),
    names_to  = "Estimator",
    values_to = "Value"
  )
alt_colours <- rep(c("#444441", "#B4B2A9"), length.out = length(populations))
names(alt_colours) <- populations

p_bar <- ggplot(diversity_long, aes(x = Population, y = Value, fill = Population)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = alt_colours) +
  facet_wrap(~ Estimator, scales = "free_y", ncol = 1) +
  labs(
    title = "Haplotype diversity estimators per population",
    x     = NULL,
    y     = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

print(p_bar)
ggsave("population_diversity_barplots.png", p_bar, width = 8, height = 13, dpi=300)


# ── 4. Haplotype diversity - pi correlation ───────────────────────────────────


