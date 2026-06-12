rm(list=ls())

# load libraries
invisible(lapply(c("here", "pegas", "tidyverse", "raster",
                   "RColorBrewer", "ggpubr", "vegan", "robust", "ggVennDiagram", "cowplot",
                   "corrplot", "UpSetR"), library, character.only = TRUE))
# -------------------------------------------------------------
# 0. Load input data and fit the environmental model
# -------------------------------------------------------------
load("rda_input_data.RData")

#RDAenv <- rda(geno ~ area + agriculture_200 + urban_200 + transparency + conductivity + pH + N + P 
#                + chl_a + float + emerge + fish 
#                + Condition(pi + Ne_LDgenome + mit_diversity + latitude + longitude), variables)

#save(RDAenv, file = "RDAenv_model.RData")
load("RDAenv_model.RData")

# -------------------------------------------------------------
# 1. Screeplot & Biplot without outliers
# -------------------------------------------------------------

screeplot(RDAenv, main="Eigenvalues of constrained axes")

# Extract scores
TAB_loci <- as.data.frame(scores(RDAenv, display = "species", scaling = 2))
TAB_var  <- as.data.frame(scores(RDAenv, display = "bp",      scaling = 2))

# Plot
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
#p_biplot
ggsave("biplot_RDAenv_RDA12.png", p_biplot, height=6, width=10, dpi=200)

# -------------------------------------------------------------
# 2, 3D biplot
# -------------------------------------------------------------
library(plotly)

# Extract scores for 3 axes
TAB_loci <- as.data.frame(scores(RDAenv, display = "species", scaling = 2, choices = 1:3))
TAB_var  <- as.data.frame(scores(RDAenv, display = "bp",      scaling = 2, choices = 1:3))

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
htmlwidgets::saveWidget(p, "biplot_3d_RDAenv.html", selfcontained = TRUE)

# -------------------------------------------------------------
# 3. Biplot with outliers
# tbc
# -------------------------------------------------------------