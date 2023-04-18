############################################################
# Description: A Script for DMP visualisation
############################################################

library("methylKit")
library("GenomicRanges")
library("genomation")
library("ggplot2")

message("Loading Annotation")
gff <- gffToGRanges("../03_DMP/Danio_rerio.GRCz11.98.gtf.gz")

message("Loading DMP.")
load("../03_DMP/Data/myDiff.rda")
myDiff15p <- getData(getMethylDiff(myDiff, difference = 15, qvalue = 0.05))

message("Draw Gene Feature Barplot")
source("./GeneFeatureBP.R")
GeneFeatureBP(myDiff15p, gff)

message("Draw Volcano Plot")
source("./VolcanoPlot.R")
VolcanoPlot(myDiff, gff)

