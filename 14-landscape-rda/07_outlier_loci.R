rm(list=ls())

# load libraries
invisible(lapply(c("here", "pegas", "tidyverse", "raster",
                   "RColorBrewer", "ggpubr", "vegan", "robust", "ggVennDiagram", "cowplot",
                   "corrplot", "UpSetR"), library, character.only = TRUE))
library(ggrepel)
library(patchwork)
library(tidyverse)
source('RDA_Functions.R')

# -------------------------------------------------------------
# 0. Load input data and fit the environmental models
# -------------------------------------------------------------
load("rda_input_data.RData")

RDAagri <- rda(geno ~ agriculture_200 + 
                 Condition(latitude + longitude + area + urban_200 + transparency 
                           + conductivity + pH + N + P + chl_a + float + emerge + fish 
                           + pi + Ne_LDgenome + mit_diversity), variables)

RDAurban <- rda(geno ~ urban_200 + 
                 Condition(latitude + longitude + area + agriculture_200 + transparency 
                           + conductivity + pH + N + P + chl_a + float + emerge + fish 
                           + pi + Ne_LDgenome + mit_diversity), variables)
summary(RDAurban)

save(RDAagri, file = "RDAagri_model.RData")
save(RDAurban, file = "RDAurban_model.RData")
load("RDAagri_model.RData")
load("RDAurban_model.RData")

### side quest begin ### 
# ======================
#to check the effect of agriculture in unconstrained vs constrained rda
RDAagrionly <- rda(geno ~ agriculture_200, variables)
summary(RDAagrionly)
summary(RDAagri)
# Results:
# Total variance = 6739.9
# agrionly = 168.6
# agri = 87.3
# --> some of the variance explained is conditioned with the other variables
### side quest end ###

# -------------------------------------------------------------
#  Preparation and manhattan plot of outlier SNPs
#   - updated rdadapt to include K=1
#   - But λ correction is not appropriate for low variance explained with K=1
#     it is sensitive to whether the axis carries enough signal for the median to meaningfully represent the null. 
#     With a weak environmental gradient or heavily conditioned model like here, it breaks down.
#  For just agriculture and urbanism I am just using the traditional chi2(1)-based p-values with the bonferroni threshold
#   - It is the most conservative and least-assuming approach --> defensible
# -------------------------------------------------------------

# =============================================================
#                       AGRICULTURE
# =============================================================

# Screeplot to explore the importance of the different RDA axes.
# Here we have only one...
screeplot(RDAagri, main="Eigenvalues of constrained axes")

# diagnostics of p-value overinflation:
# =====================
# test that my λ values are ok
zscores <- RDAagri$CCA$v[, 1]
resscale <- scale(zscores)
resmaha <- resscale^2
lambda <- median(resmaha) / qchisq(0.5, df = 1)
cat("Lambda:", lambda, "\n")  # Should be much higher, so now it overinflates all p-values

# Diagnose this
summary(eigenvals(RDAagri))   # What % variance does CCA1 explain? -> 0.0176
hist(zscores, breaks = 50)     # Are loadings extremely tight around 0? -> yes
hist(resmaha, breaks = 50)     # Are most values near 0? -> yes

# result: the p-values are tighlty around 0 --> no normality (leptokurtic distribution) --> λ is way below the useful range 
# end of diagnostics
# =====================

# Probably for loci of many environmental gradients it will be useful:
# Running the adapted function with K = 1 = number of RDA axes to keep
# rdadapt_agri<-rdadapt(RDAagri, 1)

# Instead I just use the raw chi²(1) p-values (skip λ correction)
zscores      <- RDAagri$CCA$v[, 1]
resscale     <- scale(zscores)
resmaha      <- resscale^2
pvalues_raw  <- pchisq(resmaha, df = 1, lower.tail = FALSE)

# Bonferroni threshold
thres_env <- 0.05 / length(pvalues_raw)

# Outlier identification
snp_info <- do.call(rbind, strsplit(colnames(geno), "_"))

outliers <- data.frame(
  Loci     = colnames(geno)[which(pvalues_raw < thres_env)],
  p.value  = pvalues_raw[which(pvalues_raw < thres_env)],
  scaffold = as.integer(snp_info[which(pvalues_raw < thres_env), 2]),
  position = as.integer(snp_info[which(pvalues_raw < thres_env), 3])
)

outliers_rdadapt_agri <- as.character(outliers$Loci[!duplicated(outliers$scaffold)])

# Summary statistics
total_loci   <- ncol(geno)
outlier_count <- nrow(outliers)
outlier_pct  <- round((outlier_count / total_loci) * 100, 2)

cat("Total loci:    ", total_loci, "\n")
cat("Outlier loci:  ", outlier_count, "\n")
cat("Percentage:    ", outlier_pct, "%\n")

# Manhattan plot
Outliers <- factor(ifelse(colnames(geno) %in% outliers$Loci, "Outlier", "Neutral"),
                   levels = c("Neutral", "Outlier"))

TAB_manhatan_agri <- data.frame(
  scaffold = as.integer(snp_info[, 2]),
  pos      = as.integer(snp_info[, 3]) / 1e6,
  pvalues  = pvalues_raw,
  Outliers = Outliers
)

p_agri <- ggplot(data = TAB_manhatan_agri) +
  geom_point(aes(x = pos, y = -log10(pvalues), col = Outliers), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF")) +
  labs(title = "Manhattan plot - Agriculture",
       x = "Position (Mbp)", y = "-log10(p.values)") +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed",
             color = gray(.80), linewidth = 0.6) +
  facet_wrap(~scaffold, nrow = 3, scales = "free_x") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position = "right", legend.background = element_blank(),
        panel.grid = element_blank(), legend.box.background = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = rel(.8)),
        strip.text = element_text(size = 11))


p_agri

ggsave("plots/outliers_agriculture.png", width=10, height=8, p_agri, dpi=300)

# Outlier table
outliers_table <- data.frame(
  Loci     = outliers$Loci,
  Scaffold = outliers$scaffold,
  Position = outliers$position,
  p.value  = formatC(outliers$p.value, format = "e", digits = 3)
) %>%
  arrange(Scaffold, Position)

write.csv(outliers_table,
          "results/outliers_agriculture.csv",
          row.names = FALSE)

# Site scores plot
sites_df <- data.frame(
  score      = RDAagri$CCA$u[, 1],
  agri       = variables$agriculture_200,
  population = variables$population
)

sites_pop <- sites_df %>%
  group_by(population) %>%
  summarise(score = mean(score),
            agri  = mean(agri))

p_sites_agri <- ggplot(sites_pop, aes(x = agri, y = score, label = population)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3, family = "Times",
                  box.padding = 0.4,
                  max.overlaps = Inf) +
  geom_smooth(method = "lm", se = FALSE, color = "#185FA5") +
  labs(title = "Site scores vs agriculture gradient",
       x = "Agriculture (200m)", y = "RDA1 site score") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("plots/sites_agriculture.png", width=10, height=8, p_sites_agri, dpi=300)


# =============================================================
#                        URBANISM
# =============================================================

# No need for rdadapt since K-1
# rdadapt_urban<-rdadapt(RDAurban, 1)

# Instead I just use the raw chi²(1) p-values (skip λ correction)
zscores      <- RDAurban$CCA$v[, 1]
resscale     <- scale(zscores)
resmaha      <- resscale^2
pvalues_raw  <- pchisq(resmaha, df = 1, lower.tail = FALSE)

# Bonferroni threshold
thres_env <- 0.05 / length(pvalues_raw)

# Outlier identification
snp_info <- do.call(rbind, strsplit(colnames(geno), "_"))

outliers <- data.frame(
  Loci     = colnames(geno)[which(pvalues_raw < thres_env)],
  p.value  = pvalues_raw[which(pvalues_raw < thres_env)],
  scaffold = as.integer(snp_info[which(pvalues_raw < thres_env), 2]),
  position = as.integer(snp_info[which(pvalues_raw < thres_env), 3])
)

outliers_rdadapt_urban <- as.character(outliers$Loci[!duplicated(outliers$scaffold)])

# Manhattan plot
Outliers <- factor(ifelse(colnames(geno) %in% outliers$Loci, "Outlier", "Neutral"),
                   levels = c("Neutral", "Outlier"))

TAB_manhatan_urban <- data.frame(
  scaffold = as.integer(snp_info[, 2]),
  pos      = as.integer(snp_info[, 3]) / 1e6,
  pvalues  = pvalues_raw,
  Outliers = Outliers
)

p_urban <- ggplot(data = TAB_manhatan_urban) +
  geom_point(aes(x = pos, y = -log10(pvalues), col = Outliers), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF")) +
  labs(title = "Manhattan plot - Urbanism",
       x = "Position (Mbp)", y = "-log10(p.values)") +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed",
             color = gray(.80), linewidth = 0.6) +
  facet_wrap(~scaffold, nrow = 3, scales = "free_x") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position = "right", legend.background = element_blank(),
        panel.grid = element_blank(), legend.box.background = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = rel(.8)),
        strip.text = element_text(size = 11))
p_urban

ggsave("plots/outliers_urbanism.png", width=10, height=8, p_urban, dpi=300)

# Summary statistics
total_loci   <- ncol(geno)
outlier_count <- nrow(outliers)
outlier_pct  <- round((outlier_count / total_loci) * 100, 2)

cat("Total loci:    ", total_loci, "\n")
cat("Outlier loci:  ", outlier_count, "\n")
cat("Percentage:    ", outlier_pct, "%\n")

# Outlier table
outliers_table <- data.frame(
  Loci     = outliers$Loci,
  Scaffold = outliers$scaffold,
  Position = outliers$position,
  p.value  = formatC(outliers$p.value, format = "e", digits = 3)
) %>%
  arrange(Scaffold, Position)

write.csv(outliers_table,
          "results/outliers_urbanism.csv",
          row.names = FALSE)

# Site scores plot
sites_df <- data.frame(
  score      = RDAurban$CCA$u[, 1],
  urban       = variables$urban_200,
  population = variables$population
)

sites_pop <- sites_df %>%
  group_by(population) %>%
  summarise(score = mean(score),
            urban  = mean(urban))

p_sites_urban <- ggplot(sites_pop, aes(x = urban, y = score, label = population)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3, family = "Times",
                  box.padding = 0.4,
                  max.overlaps = Inf) +
  geom_smooth(method = "lm", se = FALSE, color = "#185FA5") +
  labs(title = "Site scores vs urbanism gradient",
       x = "Urbanism (200m)", y = "RDA1 site score") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("plots/sites_urbanism.png", width=10, height=8, p_sites_urban, dpi=300)

# =============================================================
#                    ZOOM: SCAFFOLD 5, 1.1-1.3 Mbp
# =============================================================

zoom_range <- c(1.15, 1.25)
gene_pos1 <- 1200850 / 1e6
gene_pos2 <- 1201004 / 1e6

# Filter to scaffold 5 and zoom region
zoom_agri <- TAB_manhatan_agri[TAB_manhatan_agri$scaffold == 5 &
                                 TAB_manhatan_agri$pos >= zoom_range[1] &
                                 TAB_manhatan_agri$pos <= zoom_range[2], ]

zoom_urban <- TAB_manhatan_urban[TAB_manhatan_urban$scaffold == 5 &
                                   TAB_manhatan_urban$pos >= zoom_range[1] &
                                   TAB_manhatan_urban$pos <= zoom_range[2], ]

# Add population label for faceting
zoom_agri$population  <- "Agriculture"
zoom_urban$population <- "Urbanism"

zoom_combined <- rbind(zoom_agri, zoom_urban)

# Shared y-axis limits
y_max <- max(-log10(zoom_combined$pvalues), na.rm = TRUE)

p_zoom <- ggplot(data = zoom_combined) +
  geom_point(aes(x = pos, y = -log10(pvalues), col = Outliers), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF")) +
  geom_vline(xintercept = gene_pos1, color = "red", linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = gene_pos2, color = "red", linetype = "dotted", linewidth = 0.8) +
  labs(title = "Scaffold 5 — 1.1 to 1.3 Mbp",
       x = "Position (Mbp)", y = "-log10(p.values)") +
  facet_wrap(~population, nrow = 2) +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position = "right", legend.background = element_blank(),
        panel.grid = element_blank(), legend.box.background = element_blank(),
        plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold"))
p_zoom
ggsave("plots/zoom_scaffold5_1.1_1.3Mb.png", p_zoom, width = 6, height = 8, dpi = 300)

