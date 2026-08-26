rm(list = ls())

# Load libraries
library(ggplot2)

# File locations
stats_file <- "summaries/batch1234_Mt.raw.summary.txt"
plot <- "plots/05_mitochondrial/depth_distribution.png"

# Read all lines
lines <- readLines(stats_file)

# Extract DP lines
dp_lines <- grep("^DP", lines, value = TRUE)

# Convert to data frame
dp <- read.table(
  text = paste(dp_lines, collapse = "\n"),
  header = FALSE,
  stringsAsFactors = FALSE
)

# Check columns
head(dp)

# bcftools DP format is usually:
# V1=DP, V2=bin, V3=number of sites, V4=fraction
colnames(dp) <- c("DP", "tag", "depth", "sites", "fraction")

library(ggplot2)

# Make sure depth and sites are numeric
dp$depth <- as.numeric(as.character(dp$depth))
dp$sites <- as.numeric(as.character(dp$sites))

# Plot DP distribution, max read depth = 100, more x-axis detail
p <- ggplot(dp, aes(x = depth, y = sites, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  theme_bw() +
  labs(
    title = "Depth (DP) Distribution",
    x = "Read Depth (DP)",
    y = "Number of Sites"
  )
ggsave(plot, p)
