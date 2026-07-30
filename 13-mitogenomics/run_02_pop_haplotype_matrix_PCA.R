# =============================================================================
#   Produce mitohaplotype-population table
# =============================================================================

library(stringr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(tidyr)
library(dplyr)

# ── 1. Read the file ─────────────────────────────────────────────────────────
file_path <- "nexus_files/4. DnaSP_pop_haplotypes.nex"  # <-- change to your actual file path
input <- paste(readLines(file_path), collapse = "\n")

# ── 2. Extract population names ──────────────────────────────────────────────
trait_line  <- str_extract(input, "TraitLabels[^\n;]+")
populations <- str_split(str_trim(str_remove(trait_line, "TraitLabels")), "\\s+")[[1]]

# ── 3. Extract matrix lines ──────────────────────────────────────────────────
matrix_block <- str_extract(input, "(?<=Matrix\n)(.|\n)+?(?=\n;)")
matrix_lines <- str_trim(str_split(matrix_block, "\n")[[1]])
matrix_lines <- matrix_lines[nzchar(matrix_lines)]

# ── 4. Build the data frame ──────────────────────────────────────────────────
parsed <- lapply(matrix_lines, function(line) {
  parts  <- str_split(line, "\\s+", n = 2)[[1]]
  hap    <- parts[1]
  counts <- as.integer(str_split(parts[2], ",")[[1]])
  c(Haplotype = hap, setNames(counts, populations))
})

df <- do.call(rbind, lapply(parsed, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
df[populations] <- lapply(df[populations], as.integer)

# ── 5. Preview & export ──────────────────────────────────────────────────────
print(df)

write.csv(df, "pop_mitohaplotypes_matrix.csv", row.names = FALSE)

# =============================================================================
#   PCA on the haplotype x population matrix
# =============================================================================

# ── 6. Decide orientation ─────────────────────────────────────────────────────

# Extract numeric matrix (haplotypes x populations)
mat <- as.matrix(df[, populations])
rownames(mat) <- df$Haplotype

# Remove haplotypes with zero counts across all populations (e.g. Hap_40)
mat <- mat[rowSums(mat) > 0, ]

# ── 7. Normalize — convert counts to proportions within each population ───────
# This removes the effect of unequal sampling sizes across populations
mat_prop <- sweep(mat, 2, colSums(mat), FUN = "/")
mat_prop[is.nan(mat_prop)] <- 0   # fix NaN from empty populations if any

# ── 8. PCA  ────────────────────────────────────────────────────────────────────
pca_pop <- prcomp(t(mat_prop), center = TRUE, scale. = TRUE)

# Variance explained
var_pop <- round(summary(pca_pop)$importance[2,] * 100, 1)
# Score data frame for plotting
scores_pop <- as.data.frame(pca_pop$x)
scores_pop$Population <- rownames(scores_pop)

# ── 9. Scree plot — variance explained per PC ─────────────────────────────────
scree_pop <- data.frame(
  PC  = paste0("PC", 1:length(pca_pop$sdev)),
  Var = round(summary(pca_pop)$importance[2, ] * 100, 1)
)
scree_pop$PC <- factor(scree_pop$PC, levels = scree_pop$PC)

p_scree <- ggplot(scree_pop[1:10, ], aes(x = PC, y = Var)) +
  geom_col(fill = "#185FA5", width = 0.6) +
  geom_text(aes(label = paste0(Var, "%")), vjust = -0.5, size = 3) +
  labs(title = "Scree plot — populations PCA", y = "Variance explained (%)", x = "") +
  theme_classic()

print(p_scree)
ggsave("pca_mitohaplotypes_scree.png", p_scree, width = 7, height = 4, dpi=300)

# ── 10. Plot PCA, 2 PCs ───────────────────────────────────────────────────────
p_pop <- ggplot(scores_pop, aes(x = PC1, y = PC2, colour = Group, label = Population)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    title = "PCA mitohaplotypes",
    x = paste0("PC1 (", var_pop[1], "%)"),
    y = paste0("PC2 (", var_pop[2], "%)")
  ) +
  theme_classic() +
  theme(legend.position = "right")

print(p_pop)
ggsave("pca_mitohaplotypes_PC12.png", p_pop, width = 8, height = 6, dpi=300)

# ── 11. 3D PCA plot ───────────────────────────────────────────────────────────
fig_pop <- plot_ly(
  scores_pop,
  x = ~PC1, y = ~PC2, z = ~PC3,
  text  = ~Population,
  hovertemplate = paste(
    "<b>%{text}</b><br>",
    "PC1: %{x:.2f}<br>",
    "PC2: %{y:.2f}<br>",
    "PC3: %{z:.2f}<extra></extra>"
  ),
  type   = "scatter3d",
  mode   = "markers+text",
  marker = list(size = 6, color = "#185FA5")
) |>
  layout(
    title = "PCA mitohaplotypes — PC1 / PC2 / PC3",
    scene = list(
      xaxis = list(title = paste0("PC1 (", var_pop[1], "%)")),
      yaxis = list(title = paste0("PC2 (", var_pop[2], "%)")),
      zaxis = list(title = paste0("PC3 (", var_pop[3], "%)"))
    )
  )
print(fig_pop)
htmlwidgets::saveWidget(fig_pop, "pca_mithaplotypes_3D.html", selfcontained = TRUE)

# ── 12. Geography correlation ─────────────────────────────────────────────────
coord <- read.csv("metadata_files/pond_coordinates.csv")
coord

scores_geo <- scores_pop |>
  select(Population, PC1, PC2, PC3, PC4, PC5) |>
  left_join(coord, by = c("Population" = "Site"))

pcs <- paste0("PC", 1:5)

cor_results <- do.call(rbind, lapply(pcs, function(pc) {
  do.call(rbind, lapply(c("Latitude", "Longitude"), function(geo) {
    test <- cor.test(scores_geo[[pc]], scores_geo[[geo]], method = "spearman")
    data.frame(
      PC        = pc,
      Geography = geo,
      r         = round(test$estimate, 3),
      p_value   = round(test$p.value, 4),
      sig       = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01  ~ "**",
        test$p.value < 0.05  ~ "*",
        test$p.value < 0.1   ~ ".",
        TRUE                 ~ "ns"
      )
    )
  }))
}))

rownames(cor_results) <- NULL
print(cor_results)
write.csv(cor_results, "pca_geo_correlations.csv", row.names = FALSE)
