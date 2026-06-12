rm(list=ls())

# load libraries
invisible(lapply(c("here", "pegas", "tidyverse", "raster",
                   "RColorBrewer", "ggpubr", "vegan", "robust", "ggVennDiagram", "cowplot",
                   "corrplot", "UpSetR"), library, character.only = TRUE))

# -------------------------------------------------------------
# 0. Adapt .012 out files from vcftools -012 command
# -------------------------------------------------------------
geno <- read.table("data/batch1234_DmagnaLRV01_merged1_10.filtered.012.012", row.names = 1)
indv <- read.table("data/batch1234_DmagnaLRV01_merged1_10.filtered.012.012.indv")
pos  <- read.table("data/batch1234_DmagnaLRV01_merged1_10.filtered.012.012.pos")
variables <- read.csv("data/variables_samples.csv")
variables$fish<-as.factor(variables$fish)
popmap <- read.table("data/population_map.txt", sep = "", header=T)

# Add line and population names in the variables object
variables <- cbind(
  X          = variables$X,
  sample     = popmap$sample,
  population = popmap$population,
  variables[, -1]           # all columns except X
)

# Set row and column names of geno
rownames(geno) <- indv$V1
colnames(geno) <- paste(pos$V1, pos$V2, sep = "_")

save(geno, variables, file = "rda_input_data.RData")
load("rda_input_data.RData")
# -------------------------------------------------------------
# 1. Variable correlation
# -------------------------------------------------------------
library(ggcorrplot)

# Compute correlation matrix
cor_mat <- cor(variables[,(1:20)], use = "pairwise.complete.obs")
# With p-values
p_mat <- cor_pmat(variables[,(1:20)])
# plot iy & save
p <- ggcorrplot(cor_mat,
           p.mat     = p_mat,
           type      = "lower",
           lab       = TRUE,        # show correlation values
           lab_size  = 3,
           sig.level = 0.05,        # cross out non-significant
           insig     = "blank")     # or "pch" to show X
ggsave("variable_correlations.png", p, width = 9, height = 9, dpi=350)

# -------------------------------------------------------------
# 2. Ecological variable selection - Forward selection
# -------------------------------------------------------------
# Notes: 
#   - there is significant, increased correlation between longitude 
#       and agricutlure and conductivity (because of salinity)
#   - phyco~chl_a = 0.86
#   - pesticides~float = 0.83
#   - agriculture~trees = -0.73
#   - N~P = 0.63
#   - submerge~emerge = 0.69
# Solutions:
#   - I take chlor_a over phyco (they both show productivity, so I took randomly since they have similar collinearities with the rest)
#   - I take float over pesticides, as pesticides are known to be very quickly degradable (hence unreliable?) and it makes sense that the floating plants could interfere with the degradation levels of pesticides and/or water sampling.  
#   - I take agriculture over trees, as it was more precisely calculated with satelite data
#   ...
#   - I take emerge over submerge, as it is more easily identifiable (=more reliable?) plus no much difference in other interactions
# -------------------------------------------------------------

# The NULL model
RDA0 <- vegan::rda(geno ~ 1,  variables) 

# The FULL model
RDAfull <- rda(geno ~ area + agriculture_200 + urban_200 
               + transparency + conductivity + pH + N + P + chl_a
               + float + emerge + fish, variables)

# Run ordiR2step in batch 
# (running bash script that runs an R script)
system("sbatch run_ordiR2step.sh")
system("squeue -M genius")

# Load the ordiR2step results
load("var.sel.RData")
sel.anova <- read.csv("ordiR2step_anova.csv")

# -------------------------------------------------------------
### Results ###
# -------------------------------------------------------------
# all included variables in the RDAfull (after removing the highly correlated ones as mentioned)
# gave significance p=0.002 which was the lowest threshold with permutation of 1000 (?)

# -------------------------------------------------------------
# 3. Demography variable selection - Forward selection
# -------------------------------------------------------------
RDA0 <- vegan::rda(geno ~ 1,  variables) 

RDAfull <- rda(geno ~ pi + Ne_LDgenome + mit_diversity, variables)


# running forward selection of variables with ordiR2step
# The Rscript from before...
sink("ordiR2step_demo_output.txt")
var.sel <- ordiR2step(RDA0, RDAfull, Pin = 0.05, R2permutations = 1000)
sink()

capture.output(var.sel$call, file = "ordiR2step_demo_output.txt", append = TRUE)
capture.output(var.sel$anova, file = "ordiR2step_demo_output.txt", append = TRUE)

write.csv(var.sel$anova, "ordiR2step_demo_anova.csv")
