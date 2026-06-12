rm(list=ls())

# load libraries
invisible(lapply(c("here", "pegas", "tidyverse", "raster",
                   "RColorBrewer", "ggpubr", "vegan", "robust", "ggVennDiagram", "cowplot",
                   "corrplot", "UpSetR"), library, character.only = TRUE))
# -------------------------------------------------------------
# 0. Load input data and fit the environmental model
# -------------------------------------------------------------
load("rda_input_data.RData")

RDAgeog <- rda(geno ~ latitude + longitude + 
                 Condition(area + agriculture_200 + urban_200 + transparency 
                           + conductivity + pH + N + P + chl_a + float + emerge + fish 
                           + pi + Ne_LDgenome + mit_diversity), variables)

save(RDAgeog, file = "RDAgeog_model.RData")
load("RDAgeog_model.RData")

# -------------------------------------------------------------
# 1. Screeplot & Biplot of SNPs without outliers
# -------------------------------------------------------------

screeplot(RDAgeog, main="Eigenvalues of constrained axes")

# Extract scores
TAB_loci <- as.data.frame(scores(RDAgeog, display = "species", scaling = 2))
TAB_var  <- as.data.frame(scores(RDAgeog, display = "bp",      scaling = 2))

# plot
p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_point(data = TAB_loci, aes(x = RDA1 * 20, y = RDA2 * 20),
             colour = "gray70", size = 1.4) +
  geom_segment(data = TAB_var, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               colour = "red", size = 0.3, linetype = 1,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x = 1.3 * RDA1, y = 1.3 * RDA2, label = rownames(TAB_var)),
            size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        strip.text = element_text(size = 20))
p_biplot
ggsave("biplot_RDAgeog_RDA12.png", p_biplot, height=6, width=10, dpi=300)

# -------------------------------------------------------------
# 2. 3d biplot of SNPs
# -------------------------------------------------------------
library(plotly)

# Extract scores for 3 axes
TAB_loci <- as.data.frame(scores(RDAgeog, display = "species", scaling = 2, choices = 1:3))
TAB_var  <- as.data.frame(scores(RDAgeog, display = "bp",      scaling = 2, choices = 1:3))

# Scale factors
loci_scale  <- 20   # loci points
arrow_scale <- 8    # increase this to make arrows longer
label_scale <- 1.3  # label offset from arrow tip

TAB_loci_scaled <- TAB_loci * loci_scale
TAB_var_scaled  <- TAB_var  * arrow_scale

# Base plot: loci points
p <- plot_ly() %>%
  add_trace(
    data = TAB_loci_scaled,
    x = ~RDA1, y = ~RDA2, z = ~RDA3,
    type = "scatter3d", mode = "markers",
    marker = list(size = 2, color = "gray70"),
    name = "Loci",
    showlegend = TRUE
  )

# Add one trace per arrow + label
for (i in 1:nrow(TAB_var_scaled)) {
  p <- p %>%
    add_trace(
      x = c(0, TAB_var_scaled$RDA1[i]),
      y = c(0, TAB_var_scaled$RDA2[i]),
      z = c(0, TAB_var_scaled$RDA3[i]),
      type = "scatter3d", mode = "lines",
      line = list(color = "#e31a1c", width = 5),
      name = rownames(TAB_var)[i],
      showlegend = FALSE
    ) %>%
    add_trace(
      x = TAB_var_scaled$RDA1[i] * label_scale,
      y = TAB_var_scaled$RDA2[i] * label_scale,
      z = TAB_var_scaled$RDA3[i] * label_scale,
      type = "scatter3d", mode = "text",
      text = rownames(TAB_var)[i],
      textfont = list(size = 14, family = "Times", color = "#e31a1c"),
      showlegend = FALSE
    )
}

# Layout
p <- p %>%
  layout(
    scene = list(
      xaxis = list(title = "RDA 1", zeroline = TRUE),
      yaxis = list(title = "RDA 2", zeroline = TRUE),
      zaxis = list(title = "RDA 3", zeroline = TRUE)
    )
  )

p
htmlwidgets::saveWidget(p, "biplot_3d_RDAgeog.html", selfcontained = TRUE)


# -------------------------------------------------------------
# 3. Biplot of samples
# -------------------------------------------------------------

# Extract scores
TAB_sites <- as.data.frame(scores(RDAgeog, display = "sites", scaling = 2))
TAB_var   <- as.data.frame(scores(RDAgeog, display = "bp",    scaling = 2))

# To colour by population
TAB_sites$population <- variables$population
pop_longitudes <- aggregate(longitude ~ population, data = variables, FUN = mean)
pop_longitudes$lon_rank <- rank(pop_longitudes$longitude)  # added
TAB_sites <- merge(TAB_sites, pop_longitudes, by = "population")

colours <- c(
  "#8B00FF",  # 1  ultraviolet (violet)
  "#7400FF",  # 2  violet-purple
  "#5500FF",  # 3  blue-violet
  "#3300FF",  # 4  deep blue
  "#0000FF",  # 5  blue
  "#0044FF",  # 6  blue-azure
  "#0088FF",  # 7  azure
  "#00CCFF",  # 8  cyan-blue
  "#00FFFF",  # 9  cyan
  "#00FFAA",  # 10 cyan-green
  "#00FF55",  # 11 green-cyan
  "#00FF00",  # 12 green
  "#AAFF00",  # 13 yellow-green
  "#FFFF00",  # 14 yellow
  "#FFCC00",  # 15 yellow-orange
  "#FF8800",  # 16 orange
  "#FF4400",  # 17 orange-red
  "#FF0000",  # 18 red
  "#CC0000",  # 19 deep red
  "#880000"   # 20 infrared (dark red)
)

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_point(data = TAB_sites, aes(x = RDA1, y = RDA2, colour = lon_rank),
             size = 1.4) +
  scale_colour_gradientn(colours = colours, name = NULL,
                         breaks = 1:20,                                      
                         labels = pop_longitudes$population[order(pop_longitudes$longitude)]) +
  geom_segment(data = TAB_var, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               colour = "black", size = 0.3, linetype = 1,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x = 1.3 * RDA1, y = 1.3 * RDA2, label = rownames(TAB_var)),
            size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        strip.text = element_text(size = 20))

p_biplot
ggsave("biplot_RDAdemo_RDA12_sites.png", p_biplot, height = 6, width = 10, dpi = 100)


