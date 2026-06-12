rm(list=ls())

# load libraries
invisible(lapply(c("here", "pegas", "tidyverse", "raster",
                   "RColorBrewer", "ggpubr", "vegan", "robust", "ggVennDiagram", "cowplot",
                   "corrplot", "UpSetR"), library, character.only = TRUE))
library(eulerr)
library(vegan)
library(grid)
# -------------------------------------------------------------
# 0. Load input data and fit the full model
# -------------------------------------------------------------
load("rda_input_data.RData")

## Full model:
pRDAfull <- rda(geno ~ pi + Ne_LDgenome + mit_diversity 
                + latitude + longitude + area 
                + agriculture_200 + urban_200 + transparency
                + conductivity + pH + N + P + chl_a
                + float + emerge + fish, variables)

#anova(pRDAfull) ## Significance testing of the model
#summary(pRDAfull) ## Summary of the Full RDA model

# -------------------------------------------------------------
# 1. Fit partial models for demography, distance and ecology
# -------------------------------------------------------------
## Pure ecology model
pRDAecol <- rda(geno ~ area + agriculture_200 + urban_200 + transparency + conductivity + pH + N + P 
                + chl_a + float + emerge + fish 
                + Condition(pi + Ne_LDgenome + mit_diversity + latitude + longitude), variables)

#anova(pRDAecol)
#summary(pRDAecol)

## Pure demography 
pRDAdemo <- rda(geno ~ pi + Ne_LDgenome + mit_diversity 
                + Condition(area + agriculture_200 + urban_200 
                            + transparency + conductivity + pH + N + P 
                            + chl_a + float + emerge + fish + latitude + longitude), variables)

## Pure geographical distance  
pRDAgeog <- rda(geno ~ longitude + latitude
                + Condition(area + agriculture_200 + urban_200 
                            + transparency + conductivity + pH + N + P 
                            + chl_a + float + emerge + fish + pi 
                            + Ne_LDgenome + mit_diversity), variables)

## Ecology + demography
pRDAecol.demo <- rda(geno ~ area + agriculture_200 + urban_200 + transparency + conductivity + pH + N + P 
                + chl_a + float + emerge + fish + pi + Ne_LDgenome + mit_diversity 
                + Condition(latitude + longitude), variables)

## Geography + ecology
pRDAgeog.ecol <- rda(geno ~ area + agriculture_200 + urban_200 + transparency + conductivity + pH + N + P 
                + chl_a + float + emerge + fish + latitude + longitude
                + Condition(pi + Ne_LDgenome + mit_diversity), variables)

## Demography + geography
pRDAdemo.geog <- rda(geno ~ longitude + latitude + pi + Ne_LDgenome + mit_diversity
                     + Condition(area + agriculture_200 + urban_200 
                                 + transparency + conductivity + pH + N + P 
                                 + chl_a + float + emerge + fish), variables)

# save the models in .RData
save(pRDAfull, pRDAecol, pRDAdemo, pRDAgeog, pRDAecol.demo, pRDAgeog.ecol, pRDAdemo.geog, file = "pRDA_models.RData")
load("pRDA_models.RData")

# -------------------------------------------------------------
# 2. Statistics of the models
# -------------------------------------------------------------
models <- list(
  full      = pRDAfull,
  ecol      = pRDAecol,
  demo      = pRDAdemo,
  geog      = pRDAgeog,
  ecol.demo = pRDAecol.demo,
  geog.ecol = pRDAgeog.ecol,
  demo.geog = pRDAdemo.geog
)

full_r2adj <- RsquareAdj(models$full)$adj.r.squared
variance_table <- do.call(rbind, lapply(names(models), function(name) {
  m  <- models[[name]]
  s  <- summary(m)
  r2 <- RsquareAdj(m)
  data.frame(
    model           = name,
    total           = s$tot.chi,
    conditioned     = ifelse(is.null(s$partial.chi), NA, s$partial.chi),
    constrained     = s$constr.chi,
    unconstrained   = s$unconst.chi,
    R2              = r2$r.squared,
    R2.adj          = r2$adj.r.squared,
    R2.adj.perc     = r2$adj.r.squared * 100)
}))

write.csv(variance_table, "rda_variance_summary.csv", row.names = FALSE)


# Run anova in batch for all of them
#system("sbatch run_anova_pRDA.sh")
system("squeue -M genius")

# -------------------------------------------------------------
# 3. Plot - Upset based on pRDAs
# -------------------------------------------------------------
# Answers the questions:
# "After removing everything that demography and geography can explain, 
#  does ecology still explain genetic variation?"

variance_table$total       <- as.numeric(variance_table$total)
variance_table$constrained <- as.numeric(variance_table$constrained)

variance_table$expl.perc <- 
  (variance_table$constrained / variance_table$total) * 100

expressionInput <- c(
  `ecol&demo&geog` = variance_table$expl.perc[variance_table$model == "full"],
  ecol             = variance_table$expl.perc[variance_table$model == "ecol"],
  demo             = variance_table$expl.perc[variance_table$model == "demo"],
  geog             = variance_table$expl.perc[variance_table$model == "geog"],
  `ecol&demo`      = variance_table$expl.perc[variance_table$model == "ecol.demo"],
  `ecol&geog`      = variance_table$expl.perc[variance_table$model == "geog.ecol"],
  `demo&geog`      = variance_table$expl.perc[variance_table$model == "demo.geog"]
)

png("pRDA_upset.png", width = 2000, height = 1500, res = 300)

UpSetR::upset(
  fromExpression(expressionInput),
  order.by = "freq",
  set_size.show = FALSE,
  mainbar.y.label = "Explained variance (% of total)",
  sets            = c("ecol", "demo", "geog"),
)

dev.off()


library(ComplexUpset)
library(ggplot2)

# Convert to binary membership data frame (required by ComplexUpset)
df <- fromExpression(expressionInput)

# Define colors per set
set_colors <- c(ecol = "#1f78b4", demo = "#fdbf11", geog = "#e31a1c")

# Define colors per intersection (you control these manually)
intersection_colors <- c(
  "ecol"           = "#1f78b4",
  "demo"           = "#fdbf11",
  "geog"           = "#e31a1c",
  "ecol-demo"      = "#66CD00",  # blue + yellow = green (as in your Venn)
  "ecol-geog"      = "#f46d43",  # blue + red = orange
  "demo-geog"      = "#9e3d66",  # yellow + red = purple
  "ecol-demo-geog" = "#636363"   # all three = grey
)

upset(
  df,
  intersect = c("ecol", "demo", "geog"),
  base_annotations = list(
    "Explained variance (% of total)" = intersection_size(
      mapping = aes(fill = intersection),
      bar_number_threshold = 1
    )
  ),
  set_sizes = FALSE
) +
  scale_fill_manual(values = intersection_colors, guide = "none")



# -------------------------------------------------------------
# 4. Plot - Venn diagrams based on varpart
# -------------------------------------------------------------
# Answers the question:
# "How much variance is uniquely vs jointly explained?"

# perform varpart analysis (like the pRDA but automatically...)
ecol <- ~ area + agriculture_200 + urban_200 + transparency +
  conductivity + pH + N + P + chl_a + float + emerge + fish
demo <- ~ pi + Ne_LDgenome + mit_diversity
geog <- ~ latitude + longitude
vp <- varpart(
  geno,
  ecol,
  demo,
  geog,
  data = variables
)

# Plot simple venn diagram
plot(
  vp,
  Xnames = c("Environment", "Demography", "Geography"),
  col = cols,
  cex = 1.1,
  id.size = 1.1
)

# Plot pretty venn diagram
# Extract the adjusted R² fractions from varpart
# vp$part$fract rows: [A], [B], [C], [A∩B], [B∩C], [A∩C], [A∩B∩C], Residuals
fr <- vp$part$indfract[, 3]
vp$part$fract
# Named vector matching the euler() input names exactly
named_fr <- c(
  "Environment"                      = max(fr[1], 0),
  "Demography"                   = max(fr[2], 0),
  "Geography"                    = max(fr[3], 0),
  "Environment&Demography"           = max(fr[4], 0),
  "Environment&Geography"            = max(fr[6], 0),
  "Demography&Geography"         = max(fr[5], 0),
  "Environment&Demography&Geography" = max(fr[7], 0)
)

fit <- euler(named_fr)

# Plot & save with clean styling
png("varpart_venn.png", width = 9, height = 9, units = "in", res = 300, bg = "transparent")
plot(fit,
     fills      = list(fill = c("#1f78b4", "#fdbf11", "#e31a1c"), alpha = 0.4),
     edges      = list(col = "black", lwd = 2),
     labels     = list(font = 2, cex = 2.5),
     quantities = list(cex = 1.5,
                       labels = setNames(
                         paste0(round(named_fr * 100, 1), "%"),
                         names(named_fr)
                       ))
)
residuals <- round((1 - sum(fr[1:7])) * 100, 1)
grid.circle(
  x    = 0.70, y = 0.20,
  r    = unit(1.5, "cm"),
  gp   = gpar(fill = "grey92", col = "black", lwd = 1.2)
)
grid.text(
  label = paste0("Res.\n", residuals, "%"),
  x     = 0.70, y = 0.20,
  gp    = gpar(fontsize = 17, col = "grey30", fontface = "bold")
)
dev.off()


summary(RDAenv)

# svg save
svg("varpart_venn.svg", width = 9, height = 9, bg = "transparent")
plot(fit,
     fills      = list(fill = c("#1f78b4", "#fdbf11", "#e31a1c"), alpha = 0.4),
     edges      = list(col = "black", lwd = 2),
     labels     = list(font = 2, cex = 2.5),
     quantities = list(cex = 1.5,
                       labels = setNames(
                         paste0(round(named_fr * 100, 1), "%"),
                         names(named_fr)
                       ))
)
residuals <- round((1 - sum(fr[1:7])) * 100, 1)
grid.circle(
  x    = 0.70, y = 0.20,
  r    = unit(1.5, "cm"),
  gp   = gpar(fill = "grey92", col = "black", lwd = 1.2)
)
grid.text(
  label = paste0("Res.\n", residuals, "%"),
  x     = 0.70, y = 0.20,
  gp    = gpar(fontsize = 17, col = "grey30", fontface = "bold")
)
dev.off()



