### ============================================================
### Load libraries
### ============================================================
library(tidyverse)
library(data.table)
library(geosphere)
library(vegan)
library(ggtext)

rm(list = ls())

### ============================================================
### Settings
### ============================================================
pixy_dir <- "."
subfolder <- "pixy_stats"


### ============================================================
### Function to read DXY from PIXY
### ============================================================
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


### ============================================================
### LOAD & CLEAN DXY
### ============================================================
dxy_df <- read_pixy_stat("dxy", subfolder)

dxy_df <- dxy_df %>%
  filter(!is.na(avg_dxy), no_snps > 0)


### ============================================================
### LOAD & CLEAN DXY   (FIXED VERSION)
### ============================================================
dxy_df <- read_pixy_stat("dxy", subfolder)

# DXY output has no_sites, not no_snps
dxy_df <- dxy_df %>%
  filter(!is.na(avg_dxy), no_sites > 0)

### Global weighted DXY
global_weighted_dxy <- dxy_df %>%
  summarise(weighted_dxy = sum(avg_dxy * no_sites) / sum(no_sites))

global_weighted_dxy

### Pairwise DXY
pairwise_dxy_table <- dxy_df %>%
  group_by(pop1, pop2) %>%
  summarise(
    global_weighted_dxy = sum(avg_dxy * no_sites) / sum(no_sites),
    mean_unweighted_dxy = mean(avg_dxy),
    total_windows = n(),
    .groups = "drop"
  ) %>%
  arrange(pop1, pop2)

pairwise_dxy_table



### ============================================================
### GEOGRAPHIC DISTANCES
### ============================================================
coords <- read_csv("pond_coordinates.csv")

pairwise_distances <- coords %>%
  tidyr::crossing(coords %>% 
                    rename(Site2 = Site,
                           Latitude2 = Latitude,
                           Longitude2 = Longitude)) %>%
  filter(Site != Site2) %>%
  mutate(
    distance_m = geosphere::distHaversine(
      cbind(Longitude, Latitude),
      cbind(Longitude2, Latitude2)
    ),
    distance_km = distance_m / 1000
  ) %>%
  select(Site, Site2, distance_km) %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Site, Site2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(Site, Site2, distance_km)


### Distance matrix
dist_matrix <- pairwise_distances %>%
  pivot_wider(names_from = Site2, values_from = distance_km) %>%
  as.data.frame()

rownames(dist_matrix) <- dist_matrix$Site
dist_matrix$Site <- NULL
dist_matrix <- as.matrix(dist_matrix)
dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]


### ============================================================
### MERGE DXY WITH DISTANCE
### ============================================================
dxy_dist <- pairwise_dxy_table %>%
  rename(Site = pop1, Site2 = pop2) %>%
  inner_join(pairwise_distances, by = c("Site", "Site2"))


### ============================================================
### DXY MATRIX FOR MANTEL TEST
### ============================================================
dxy_matrix <- pairwise_dxy_table %>%
  select(pop1, pop2, global_weighted_dxy) %>%
  pivot_wider(names_from = pop2, values_from = global_weighted_dxy) %>%
  as.data.frame()

rownames(dxy_matrix) <- dxy_matrix$pop1
dxy_matrix$pop1 <- NULL

dxy_matrix <- as.matrix(dxy_matrix)
dxy_matrix[lower.tri(dxy_matrix)] <- t(dxy_matrix)[lower.tri(dxy_matrix)]


### ============================================================
### MANTEL TEST (DXY vs geographic distance)
### ============================================================
mantel_dxy <- mantel(dxy_matrix, dist_matrix, method = "pearson", permutations = 9999)
mantel_dxy

mantel_r_dxy <- round(mantel_dxy$statistic, 3)
mantel_p_dxy <- signif(mantel_dxy$signif, 3)


### ============================================================
### SCATTERPLOT WITH MANTEL ANNOTATION
### ============================================================
annot_dxy <- paste0(
  "Mantel test:\n",
  "r = ", mantel_r_dxy, "\n",
  "p = ", mantel_p_dxy
)

p_dxy <- ggplot(dxy_dist, aes(distance_km, global_weighted_dxy)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("label",
           x = Inf, y = Inf,
           label = annot_dxy,
           hjust = 1.1, vjust = 1.5,
           size = 6,
           fill = "white", color = "black") +
  theme_bw() +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise global DXY",
    title = "Isolation by Distance (DXY)"
  )

p_dxy

ggsave("Dxy_vs_Dist_Mantel.png", p_dxy, width = 8, height = 6, dpi = 300)

      