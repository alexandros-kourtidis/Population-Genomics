#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r -p <prefix> -i <info file, 2-column file with ind name and population/species name> 
#                        -k <max K value> -l <comma-separated list of populations/species in the order to be plotted>
# This R script makes barplots for K=2 and all other K values until max K (specified with -k). It labels the individuals 
# and splits them into populations or species according to the individual and population/species names in the 2-column file specified with -i.
# The order of populations/species follows the list of populations/species given with -l.
# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3
# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3
#
# Author: Joana Meier, September 2019. Adapted by Alexandros Kourtidis, May 2025

# Manually set arguments for running in RStudio
args <- list(
  prefix = "batch1234_nmaf",
  infofile = "batch1234_infofile.txt",
  minK = 13,
  maxK = 17,
  populations = "A1_AnOudin,A2_BuSN,A3_Mech,A4_Pl15YEL,A5_GenM,B1_BKN1,B2_OM2,B3_ZW,B4_OHZ,B5_DA2,C1_BlfN,C2_MO,C3_Ter1,C4_BW48630,C5_BW36962,D1_CBOO6,D2_LRV,D3_BKLE5,D4_BW62256,D5_BW22050",
  outPrefix = "batch1234_nmaf_K13-K17"
)
opt <- args

# Load optparse
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library("optparse")

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
}else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
}else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels<-read.table(opt$infofile)

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))

# Set colour palete
cluster_cols = c(
  "#F0F8FF", "#CD950C", "#BBFFFF", "#7FFFD4", "#1E90FF", "#F5F5DC",
  "#FFE4C4", "#8B8B83", "#FFF0F5", "#104E8B", "#8A2BE2", "#A52A2A",
  "#DEB887", "#5F9EA0", "#98FB98", "#D2691E", "#FF7F50"
)

# Prepare spaces to separate the populations/species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
spaces<-spaces[-length(spaces)]

# Calculate dynamic height based on number of K values
height <- 300 * (maxK - minK + 1)

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),
     width = 2000, 
     height = 1200,
     res=200)
par(mfrow=c(maxK - minK + 1,1),
    mar=c(0,1,0,0),
    oma=c(2,1,9,1),
    mgp=c(0,0.2,0),
    xaxs="i",
    cex.lab=0.8,cex.axis=0.8)

# Plot minK
bp<-barplot(t(as.matrix(tbl[[1]][order(labels$n),])), 
            col=cluster_cols[1:minK],
            xaxt="n", 
            border=NA,
            ylab=paste0("K=",minK),
            yaxt="n",
            space=spaces)
axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)

# Plot higher K values
if(maxK>minK)lapply(2:(maxK - minK + 1), function(i){
  k <- minK + i - 1
  barplot(t(as.matrix(tbl[[i]][order(labels$n),])), 
          col=cluster_cols[1:k],
          xaxt="n", border=NA,
          ylab=paste0("K=",k),
          yaxt="n",space=spaces)
})

# Add populations labels 
axis(1,
     at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
     labels=unlist(strsplit(opt$populations,",")),
     cex.axis = 0.8)

dev.off()


### Plot optimal K

# Reduce margins (especially bottom and top)
par(mar = c(10, 4, 10, 1)) # c(bottom, left, top, right)

# Set colour palete
cluster_cols = c(
  "#F0F8FF", "#CD950C", "#BBFFFF", "#7FFFD4", "#1E90FF", "#F5F5DC",
  "#FFE4C4", "#8B8B83", "#FFF0F5", "#104E8B", "#8A2BE2", "#A52A2A",
  "#DEB887", "#5F9EA0", "#98FB98", "#D2691E", "#FF7F50"
)

# Open svg device
svg("admixture_barplot_K17.svg", width = 12, height = 3)

# Plot K=17
bp<-barplot(t(as.matrix(tbl[[14]][order(labels$n),])), 
            col = cluster_cols,
            xaxt = "n", 
            border = NA,
            ylab = paste0("K=17"),
            yaxt = "n",
            space = spaces)

# Add individual labels
axis(3,at=bp,labels=labels$ind[order(labels$n)],
     las = 2,
     tick = F,
     cex.axis = 0.6)

# Add populations labels 
axis(1,
     at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
     labels=unlist(strsplit(opt$populations,",")),
     cex.axis = 0.8)

# Close svg device
dev.off()

