# plot admixture results in as pie charts on a map
# mapmixture from https://cran.r-project.org/web/packages/mapmixture/readme/README.html

# Install and load mapmixture
install.packages("mapmixture") # Run only once
library(mapmixture)
library(ggplot2)
library(sf)
library(svglite)

# Read admixture proportions
Q <- read.table("batch1234_nmaf.17.Q")

# Read infofile
info <- read.table("batch1234_infofile.txt", 
                   header = FALSE, 
                   col.names = c("Ind", "Site"))

# Combine info and admixture proportions
admixture_df <- cbind(info, Q)

# Reorder columns 1 and 2 to proper input format
admixture_df <- admixture_df[, c("Site", "Ind", setdiff(names(admixture_df), c("Site", "Ind")))]

# Rename cluster columns
colnames(admixture_df)[3:ncol(admixture_df)] <- paste0("Cluster", 1:(ncol(admixture_df) - 2))
head(admixture_df)
summary(admixture_df)

# Read in coordinates file
coordinates <- read.csv("coordinates.csv")
head(coordinates)


cluster_cols = c(
  "#F0F8FF", "#CD950C", "#BBFFFF", "#7FFFD4", "#1E90FF", "#F5F5DC",
  "#FFE4C4", "#8B8B83", "#FFF0F5", "#104E8B", "#8A2BE2", "#A52A2A",
  "#DEB887", "#5F9EA0", "#98FB98", "#D2691E", "#FF7F50"
)

# Run mapmixture
map1 <- mapmixture(admixture_df, 
                   coordinates, 
                   crs = 31370, 
                   pie_size = 0.1, 
                   cluster_cols = cluster_cols)


# Convert to sf object with WGS84
coords_sf <- st_as_sf(coordinates, coords = c("Longitude", "Latitude"), crs = 4326)

# Transform to EPSG:31370
coords_sf <- st_transform(coords_sf, crs = 31370)

# Extract X and Y for plotting
coords_sf$X <- st_coordinates(coords_sf)[, 1]
coords_sf$Y <- st_coordinates(coords_sf)[, 2]

# Add site labels above each pie chart
map1 <- map1 +
  geom_text(
    data = coords_sf,
    mapping = aes(x = X, y = Y + 10000, label = Site), # Adjust Y offset as needed
    size = 4,
    fontface = "bold"
  )

# Display the map
print(map1)

# Save the map as a png
ggplot2::ggsave("admixture_map.png", plot = map, width = 10, height = 8, dpi = 300)

# Save the map as an SVG
ggplot2::ggsave("admixture_map_labels.svg", plot = map1, width = 10, 
                height = 8, dpi = 300, device = "svg")

               