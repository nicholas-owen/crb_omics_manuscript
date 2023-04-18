############################################################
# Aurthor : Tian
# Description: An important Script, it matches the result of methylKit to gff format annotation from Ensembl
############################################################
library("methylKit")
library("genomation")
library("GenomicRanges")

# load("../3.DMP/Data/myDiff.rda")
# load("../3.DMP/Data/filteredBeta.rda")
# 
# message("Prepare DMP (cutoff 15 meth.diff and q vlaue 0.05)")
# diff <- getData(getMethylDiff(myDiff, difference=15, qvalue=0.05))
# message("!Important, keey beta matrix filtered with DMP")
# rownames(filteredBeta) <- rownames(myDiff)
# beta <- filteredBeta[rownames(diff),]
# 
# message("Loading Annotation")
# gff <- gffToGRanges("../3.DMP/Danio_rerio.GRCz11.98.gtf.gz")
# grl <- as(split(gff, gff$type), "GRangesList")
# grl[["promoter"]] <- promoters(grl[["gene"]])
# source("../4.DMPVis/ExpandGeneRange.R")
# expandedGenes <- ExpandGeneRange(grl[["gene"]])
# 
# message("To run MatchingGene function, 3 thins are required: diff, matched beta, a GRange Annotation")
# 
MatchingGene <- function(diff, beta, anno)
{

ov <- suppressWarnings(as.data.frame(findOverlaps(anno, makeGRangesFromDataFrame(diff))))
cpgnumber <- aggregate(ov$subjectHits, by = list(anno$gene_id[ov$queryHits]), function(x) length(x))
MeanMethDiff <- aggregate(ov$subjectHits, by = list(anno$gene_id[ov$queryHits]), function(x) mean(diff[x, "meth.diff"]))
MeanBeta <- aggregate(ov$subjectHits, by = list(anno$gene_id[ov$queryHits]), function(x) { if(length(x) > 1) colMeans(beta[x,]) else beta[x,]})

rownames(cpgnumber) <- cpgnumber[, 1]
rownames(MeanMethDiff) <- MeanMethDiff[, 1]
rownames(MeanBeta) <- MeanBeta[, 1]

Genes <- as.data.frame(anno[match(rownames(cpgnumber), anno$gene_id), c("gene_name", "gene_id")])
Genes$EnrichNumber <- cpgnumber[Genes$gene_id, 2]
Genes$MeanMethDiff <- MeanMethDiff[Genes$gene_id, 2]
Genes$MeanBeta <- MeanBeta[Genes$gene_id, 2]
rownames(Genes) <- Genes$gene_id

DMP <- diff
DMP$ovIndex <- 1:nrow(DMP)
tmpGene <- as.data.frame(anno[ov$queryHits, c("gene_id", "gene_name")])
tmpGene$ovIndex <- ov$subjectHits
DMP <- merge(DMP, tmpGene, by="ovIndex", all = T)

return(list(Genes=Genes, DMP=DMP))
}
