#==========
#2D and 3D visualisation of the PCA on 400 samples from 20 populations
#==========

#set my working directory
  setwd("C:/Users/u0167113/Desktop/popgen analyses/PCA/plink")

#clear the environment
rm(list = ls())

# Load necessary libraries
library("ggplot2")
library(plotly)

# Read PLINK PCA output
eigenvec_nmaf <- read.table("batch1234_nmaf.eigenvec", header = FALSE, stringsAsFactors = FALSE)
eigenval_nmaf <- scan("batch1234_nmaf.eigenval")
varprop_nmaf <- (eigenval_nmaf/sum(eigenval_nmaf))*100
varprop_nmaf

length(eigenvec_nmaf$V1)
length(eigenvec_nmaf$V2)
summary(eigenvec_nmaf)

# Extract sample names and principal components
sNames <- eigenvec_nmaf$V1
sNames <- sub("_aligned.sorted.nd.bam", "", sNames)
PCs <- eigenvec_nmaf[, -c(1, 2)] # Remove first two columns (FID and IID)

# Load population information
popInfo <- read.table("batch1234chat_metadata.txt", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
summary(popInfo)

# Filter population data to include only the samples present in the PCA in the same order as PCA sample IDs
popInfo_filtered <- popInfo[popInfo$Clonal_line_name %in% sNames, ]
popInfo_filtered <- popInfo_filtered[match(sNames, popInfo_filtered$Clonal_line_name), ]
summary(popInfo_filtered)

# Define your custom colors for each population
gradient_colors <- c("A1_AnOudin" = "#CAFF70", "A2_BuSN" = "#B3EE3A", "A3_Mech" = "#BCEE68", "A4_Pl15YEL" = "#9ACD32", "A5_GenM" = "#A2CD5A",
                   "B1_BKN1" = "#98F5FF", "B2_OM2" = "#BFEFFF", "B3_ZW" = "#B0E2FF", "B4_OHZ" = "#87CEFA", "B5_DA2" = "#A4D3EE",
                   "C1_BlfN" = "#FFF68F", "C2_MO" = "#FFEC8B", "C3_Ter1" = "#FFD700", "C4_BW48630" = "#EEC900", "C5_BW36962" = "#FFC125",
                   "D1_CBOO6" = "#FFE4E1", "D2_LRV" = "#FFB6C1", "D3_BKLE5" = "#FFAEB9", "D4_BW62256" = "#EEA2AD", "D5_BW22050" = "#CD8C95" )


# Create a data frame for plotting PCA results
plot_data <- as.data.frame(PCs)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity
plot_data$LineName <- as.character(popInfo_filtered$Clonal_line_name)
plot_data$Population <- factor(popInfo_filtered$Pond_name)
plot_data$LineNumber <- as.character(popInfo_filtered$Line_number)
plot_data$PondType <- factor(popInfo_filtered$Pond_type)
plot_data$Color <- gradient_colors[plot_data$Population]
head(plot_data)

#========  
##### Make a 2D PCA Plot with PC1 and PC2
#========

#PCA Plot Colored by Population
p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 2) +
  labs(title = "PCA Colours by Population", x = "PC1 (6.93%)", y = "PC2 (5.52%)") +
  theme_minimal() +
  theme(legend.title = element_blank())

p + scale_colour_manual(values = gradient_colors)


#========  
##### Make an interactive 3D PCA Plot with PC1, PC2, and PC3
#========

### both scatter plot and 3d meshes are interactive and with pond code annotation -

# invert the PC3 axis, optional: in order to match the orientatio of a previous snpgdsn PCA
plot_data$PC3 <- -plot_data$PC3

# Empty fig
fig <- plot_ly()

# Add volume areas for each population
for (i in seq_along(levels(plot_data$Population))) {
  pop <- levels(plot_data$Population)[i]
  pop_data <- subset(plot_data, Population == pop)
  fig <- fig %>%
    add_mesh(
      x = pop_data$PC1,
      y = pop_data$PC2,
      z = pop_data$PC3,
      opacity = 0.2,
      facecolor = rep(unique(pop_data$Color), ncol(pop_data)), 
      name = paste(pop),
      legendgroup = paste(pop),
      alphahull = 0,  # Use alphahull to create a convex hull around the points
      showscale = FALSE  
    )
}

# Add labels and title
fig <- fig %>%
  layout(
    title = "3D PCA Plot (PC1, PC2, PC3)",
    scene = list(
      xaxis = list(title = "PC1 (6.93%)"),
      yaxis = list(title = "PC2 (5.52)"),
      zaxis = list(title = "PC3 (4.60%)")
    )
  )

# Add scatter plot on top of the existing figure
# Note: The part hovering names of clonal lines 
# ( text = ~paste("Line:", LineName, "<br>Number:", LineNumber, "<br>Population:", Population) )
# only works well if the volume areas are not applied !!
fig <- fig %>%
  add_trace(
    data = plot_data, 
    x = ~PC1, y = ~PC2, z = ~PC3, 
    type = "scatter3d", 
    mode = "markers",
    marker = list(size = 4, color = plot_data$Color),
    legendgroup = ~Population,
    text = ~paste("Line:", LineName, "<br>Number:", LineNumber, "<br>Population:", Population),
    hoverinfo = "text"
  )

# Add text annotations for each population cluster
for (i in seq_along(levels(plot_data$Population))) {
  pop <- levels(plot_data$Population)[i]
  pop_data <- subset(plot_data, Population == pop)
  fig <- fig %>%
    add_text(
      x = mean(pop_data$PC1),
      y = mean(pop_data$PC2) + 0.005,
      z = mean(pop_data$PC3),
      text = substr(pop, 1, 2),
      showarrow = FALSE,
      font = list(size = 12, color = 'black')
    )
}

# Show plot
fig

