rm(list = ls())
# Load libraries
library(ggplot2)

# File locations
stats_file <- "summaries/batch1234_Mt.raw.summary.txt"
plot <- "plots/05_mitochondrial/qual_distribution.png"

# Read all lines
lines <- readLines(stats_file)
# Extract QUAL lines
qual_lines <- grep("^QUAL", lines, value = TRUE)
# Convert to data frame
qual <- read.table(
  text = paste(qual_lines, collapse = "\n"),
  header = FALSE,
  stringsAsFactors = FALSE
)
# Check columns
head(qual)

# bcftools QUAL format is usually:
# V1=QUAL, V2=id, V3=Quality, V4=number of SNPs, V5=number of transitions (1st ALT),
# V6=number of transversions (1st ALT), V7=number of indels
colnames(qual) <- c("QUAL", "tag", "quality", "snps", "ts", "tv", "indels")
library(ggplot2)
# Make sure quality and snps are numeric
qual$quality <- as.numeric(as.character(qual$quality))
qual$snps    <- as.numeric(as.character(qual$snps))
# Plot QUAL distribution
p <- ggplot(qual, aes(x = quality, y = snps, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  theme_bw() +
  labs(
    title = "Variant Quality (QUAL) Distribution",
    x = "Variant Quality (QUAL)",
    y = "Number of SNPs"
  )
ggsave(plot, p)