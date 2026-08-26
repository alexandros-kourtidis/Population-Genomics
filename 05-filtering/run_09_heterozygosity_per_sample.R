rm(list = ls())
# Load libraries
library(ggplot2)
library(stringr)
library(dplyr)

# File locations
stats_file <- "summaries/batch1234_Mt.raw.summary.txt"
plot       <- "plots/heterozygosity_per_sample.Mt_raw.png"

dir.create(dirname(plot), showWarnings = FALSE, recursive = TRUE)

# Read all lines
lines <- readLines(stats_file, warn = FALSE)

# Extract PSC lines (per-sample counts; requires bcftools stats -s -)
psc_lines <- grep("^PSC\t", lines, value = TRUE)
if (length(psc_lines) == 0) {
  stop("No PSC lines found -- was this summary file generated with `bcftools stats -s -`?")
}

# Parse PSC block
# bcftools PSC columns:
# PSC, id, sample, nRefHom, nNonRefHom, nHets, nTransitions, nTransversions,
# nIndels, avgDepth, nSingletons, nHapRef, nHapAlt, nMissing
psc <- do.call(rbind, lapply(psc_lines, function(x) {
  f <- str_split(x, "\t")[[1]]
  data.frame(
    sample     = f[3],
    nRefHom    = as.numeric(f[4]),
    nNonRefHom = as.numeric(f[5]),
    nHets      = as.numeric(f[6]),
    nMissing   = as.numeric(f[14]),
    stringsAsFactors = FALSE
  )
}))

# Heterozygosity level per sample:
# fraction of called genotypes (excluding missing) that are heterozygous
psc <- psc %>%
  mutate(
    total_called = nRefHom + nNonRefHom + nHets,
    het_frac     = nHets / total_called,
    sample       = factor(sample, levels = sort(unique(sample), decreasing = TRUE))
  )

head(psc)

# Plot heterozygosity per sample, ordered alphabetically
p <- ggplot(psc, aes(x = sample, y = het_frac)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Per-sample heterozygosity level",
    x = "Sample",
    y = "Fraction of het genotypes (nHets / total called)"
  )
ggsave(plot, p, width = 7, height = max(4, 0.15 * nrow(psc)), dpi = 100, limitsize = FALSE)
