### ============================================
### Load libraries
### ============================================
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)

rm(list = ls())

### ============================================
### Set your PIXY output directory
### ============================================
pixy_dir <- "."   # <-- change if needed

### ============================================
### Function to read all PIXY files for a stat
### ============================================
read_pixy_stat <- function(stat, subfolder = "pixy_stats") {
  stat_dir <- file.path(pixy_dir, subfolder)
  files <- list.files(
    stat_dir,
    pattern = paste0(stat, ".txt$"),
    full.names = TRUE
  )
  dt <- data.table::rbindlist(lapply(files, data.table::fread))
  dt$stat <- stat
  return(dt)
}

### ============================================
### Load all five PIXY stats
### ============================================
fst_df    <- read_pixy_stat("fst", "pixy_stats")
#dxy_df    <- read_pixy_stat("dxy", "pixy_stats")


### ============================================
### GLOBAL & PAIRWISE FST
### ============================================

# Remove windows with NA or no SNPs
fst_df <- fst_df %>% 
  filter(!is.na(avg_wc_fst), no_snps > 0)


# Global weighted FST across all population pairs and all windows
global_weighted_fst <- fst_df %>%
  summarise(
    weighted_fst = sum(avg_wc_fst * no_snps) / sum(no_snps)
  )
global_weighted_fst


# Compute global pairwise weighted FST for each pop1–pop2 pair
pairwise_fst_table <- fst_df %>%
  group_by(pop1, pop2) %>%
  summarise(
    global_weighted_fst = sum(avg_wc_fst * no_snps) / sum(no_snps),
    mean_unweighted_fst = mean(avg_wc_fst),
    total_windows = n(),
    .groups = "drop"
  ) %>%
  arrange(pop1, pop2)

pairwise_fst_table
fwrite(pairwise_fst_table, "pairwise_global_fst.csv")

# ============================================
# HEATMAP: Pairwise FST
# ============================================
library(ggplot2)
library(dplyr)
library(reshape2)

# Prepare symmetric matrix
fst_sym <- pairwise_fst_table %>%
  select(pop1, pop2, global_weighted_fst)

# Create mirrored version
fst_sym_rev <- fst_sym %>%
  rename(pop1 = pop2, pop2 = pop1)  

# Combine & remove duplicates
fst_full <- bind_rows(fst_sym, fst_sym_rev) %>%
  distinct(pop1, pop2, .keep_all = TRUE)

# Convert to matrix form
fst_wide <- fst_full %>%
  pivot_wider(names_from = pop2, values_from = global_weighted_fst) %>%
  as.data.frame()

rownames(fst_wide) <- fst_wide$pop1
fst_wide$pop1 <- NULL

# Order the axis
fst_wide <- fst_wide[sort(rownames(fst_wide)), rev(sort(colnames(fst_wide)))]

# Convert to long format for ggplot
fst_long <- melt(as.matrix(fst_wide))

# Plot heatmap
p_fst <- ggplot(fst_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = c("lightyellow", "orange", "darkred"),
    name = "FST") +
  theme_minimal() +
  labs(title = "Pairwise FST Heatmap", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank()
  )
p_fst

ggsave("pairwise_fst_heatmap.png", p_fst, width = 8, height = 6, dpi = 300)


########################
# Calculate pop distances
#########################
# Load packages
library(geosphere)   # for distHaversine
library(dplyr)
library(tidyr)
library(readr)

# 1. Read the coordinate file
coords <- read_csv("pond_coordinates.csv")

# Inspect
print(coords)

# 2. Calculate all pairwise distances (in km)
pairwise_distances <- coords %>%
  # create all combinations of site pairs
  tidyr::crossing(coords %>% rename(Site2 = Site,
                                    Latitude2 = Latitude,
                                    Longitude2 = Longitude)) %>%
  # remove identical pairs
  filter(Site != Site2) %>%
  # compute Haversine distance
  mutate(
    distance_m = geosphere::distHaversine(
      cbind(Longitude, Latitude),
      cbind(Longitude2, Latitude2)
    ),
    distance_km = distance_m / 1000
  ) %>%
  select(Site, Site2, distance_km)

# 3. Remove duplicate mirrored pairs (A–B same as B–A)
pairwise_distances <- pairwise_distances %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Site, Site2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(Site, Site2, distance_km)

# View the result
pairwise_distances

# Distance matrix
distance_matrix <- pairwise_distances %>%
  pivot_wider(names_from = Site2, values_from = distance_km)

distance_matrix

#######################################
### HEATMAP: Pairwise Geographic Distance (km)
#######################################

# Mirror the distance table
dist_sym <- pairwise_distances
dist_sym_rev <- dist_sym %>%
  rename(Site = Site2, Site2 = Site)

dist_full <- bind_rows(dist_sym, dist_sym_rev) %>%
  distinct(Site, Site2, .keep_all = TRUE)

# Wide matrix
dist_wide <- dist_full %>%
  pivot_wider(names_from = Site2, values_from = distance_km) %>%
  as.data.frame()

rownames(dist_wide) <- dist_wide$Site
dist_wide$Site <- NULL

# Order the axis
dist_wide <- dist_wide[sort(rownames(dist_wide)), rev(sort(colnames(dist_wide)))]

# Long format
dist_long <- melt(as.matrix(dist_wide))

# Plot heatmap
p_dist <- ggplot(dist_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = c("lightyellow", "orange", "darkred"),
    name = "FST") +
  theme_minimal() +
  labs(title = "Geographic Distance Heatmap", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank()
  )

ggsave("pairwise_distance_heatmap.png", p_dist, width = 8, height = 6, dpi = 300)


##############################
# Fst vs Distance plot
#############################

library(dplyr)

fst_dist <- pairwise_fst_table %>%
  rename(Site = pop1, Site2 = pop2) %>%
  inner_join(pairwise_distances, by = c("Site", "Site2"))

library(ggplot2)
summary(fst_dist)

ggplot(fst_dist, aes(x = distance_km, y = global_weighted_fst)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise global FST",
    title = "Isolation-by-distance: FST vs Geographic Distance"
  )


##############################
# Fst vs Distance plot WITH significant pairs highlighted
#############################
sig_pairs <- read.csv("significant_pairs.csv")
sig_pairs <- sig_pairs %>%
  mutate(pair_id = paste(pond1, pond2, sep = "_"))

fst_dist2 <- fst_dist %>%
  mutate(
    Site  = substr(Site, 1, 2),
    Site2 = substr(Site2, 1, 2),
    pair_id = paste(pmin(Site, Site2), pmax(Site, Site2), sep = "_")
  )

fst_dist2 <- fst_dist2 %>%
  left_join(sig_pairs %>% select(treatment, pair_id, prob_diff), 
            by = "pair_id") %>%
  mutate(is_sig = !is.na(prob_diff))   # TRUE if this pair is significant

fst_dist2 <- fst_dist2 %>%
  distinct(pair_id, .keep_all = TRUE)

fst_dist2 <- fst_dist2 %>%
  mutate(label = ifelse(is_sig, paste(Site, Site2, sep = "-"), NA))

library(ggplot2)
library(ggrepel)

p_high <- ggplot(fst_dist2, aes(x = distance_km, y = global_weighted_fst)) +
  
  # background (non-significant)
  geom_point(data = subset(fst_dist2, !is_sig),
             color = "grey70", size = 3, alpha = 0.6) +
  
  # highlighted significant points
  geom_point(data = subset(fst_dist2, is_sig),
             color = "red", size = 4, alpha = 0.9) +
  
  # labels for significant points
  ggrepel::geom_text_repel(
    data = subset(fst_dist2, is_sig),
    aes(label = label),
    size = 4,
    color = "red",
    fontface = "bold",
    max.overlaps = Inf
  ) +
  
  # regression line
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  
  theme_bw(base_size = 14) +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise global FST",
    title = "Isolation-by-distance: FST vs Geographic Distance",
    subtitle = "Red points = tukey-significant pond comparisons all treats (p < 0.05)"
  )

p_high


ggsave("fst_dist_pheno_pairwise.png", p_high, width = 8,
       height = 8,
       dpi = 300)























##############################
# Interactive Fst vs Distance plot
#############################
library(plotly)
library(htmlwidgets)

fst_dist <- fst_dist %>%
  mutate(hover_text = paste0(
    "Pair: ", Site, " – ", Site2, "\n",
    "Distance: ", round(distance_km, 2), " km\n",
    "FST: ", round(global_weighted_fst, 4)
  ))
p_int <- ggplot(fst_dist, aes(x = distance_km, y = global_weighted_fst)) +
  geom_point(aes(text = hover_text), size = 3, alpha = 0.7) +   # <-- hover text
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise global FST",
    title = "Isolation-by-distance: FST vs Geographic Distance"
  )

# Convert to interactive plotly
ggplotly(p_int, tooltip = "text")
p_int

save_html(p_int, "fst_vs_distance_interactive.html")

#######################################
# Isolation-by-distance Mantel test
#########################################
library(vegan)
library(tidyr)
library(dplyr)

fst_matrix <- pairwise_fst_table %>%
  select(pop1, pop2, global_weighted_fst) %>%
  pivot_wider(names_from = pop2, values_from = global_weighted_fst) %>%
  as.data.frame()

rownames(fst_matrix) <- fst_matrix$pop1
fst_matrix$pop1 <- NULL

fst_matrix <- as.matrix(fst_matrix)

# geographic matrix
dist_matrix <- pairwise_distances %>%
  select(Site, Site2, distance_km) %>%
  pivot_wider(names_from = Site2, values_from = distance_km) %>%
  as.data.frame()

rownames(dist_matrix) <- dist_matrix$Site
dist_matrix$Site <- NULL

dist_matrix <- as.matrix(dist_matrix)

# Symmetrize
dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]

# Symmetrize (FST is symmetric)
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]


#Mantel test

mantel_result <- mantel(fst_matrix, dist_matrix, method = "pearson", permutations = 9999)
mantel_result

################################
# Scatterplot WITH Mantel test results
################################

# Extract Mantel r and p
mantel_r <- round(mantel_result$statistic, 3)
mantel_p <- signif(mantel_result$signif, 3)

# option 1
ggplot(fst_dist, aes(x = distance_km, y = global_weighted_fst)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           label = paste0("Mantel r = ", mantel_r, "\n",
                          "p = ", mantel_p),
           size = 5) +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise global FST",
    title = "Isolation by Distance"
  )

# option 2
annot <- paste0("Mantel test:\n",
                "r = ", mantel_r, "\n",
                "* p = ", mantel_p)


p2 <- ggplot(fst_dist, aes(distance_km, global_weighted_fst)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  annotate("label",
           x = Inf, y = Inf,
           label = annot,
           hjust = 1.1, vjust = 1.5,
           size = 6,
           fill = "white", color = "black") +
  theme_bw()
p2

ggsave(filename = "Fst_vs_Dist_Mantel.png", plot = p2,
       width = 8, height = 6, dpi = 300)
