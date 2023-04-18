############################################################
# Aurthor : Tian
# Description: A Script for Integration between Methylation and RNA-seq
############################################################

library("methylKit")
library("GenomicRanges")
library("genomation")
library("ggplot2")
library("ggrepel")

message("Prepare RNA-seq Results")
DEG <- read.csv("./External/deseq2-HF-M56-annotated.diffexpr-results.csv", as.is=TRUE)
rownames(DEG) <- DEG$gene_id

message("Prepare DMP list")
load("../3.DMP/Data/myDiff.rda")
diff <- getData(myDiff)

message("Prepare Methylation Beta Value")
load("../3.DMP/Data/filteredBeta.rda")
beta <- filteredBeta

message("Loading Annotation")
gff <- gffToGRanges("../3.DMP/Danio_rerio.GRCz11.98.gtf.gz")
grl <- as(split(gff, gff$type), "GRangesList")

message("Mapping Location between Methylatio to Expression")
anno <- grl[['gene']]
message("Filter Gene Annotation with DEG list")
anno <- anno[which(anno$gene_id %in% rownames(DEG))]

Promoters <- promoters(anno, upstream=2000, downstream=0)
source("../4.DMPVis/ExpandGeneRange.R")
expandedGenes <- ExpandGeneRange(anno)

ovPromoter <- suppressWarnings(as.data.frame(findOverlaps(Promoters, as(diff, "GRanges"))))
ov <- suppressWarnings(as.data.frame(findOverlaps(expandedGenes, as(diff, "GRanges"))))

message("Form log2FC and meth.diff for Plotting")
df <- data.frame(log2FC=DEG[as.data.frame(anno)[ov$queryHits, 'gene_id'], "log2FoldChange"],
                 MethDiff=(rowMeans(beta[,4:6]) - rowMeans(beta[,1:3]))[ov$subjectHits],
                 GeneName=DEG[as.data.frame(anno)[ov$queryHits, 'gene_id'], "gene_name"]
                )

df$status <- 'No-Sig'
df$status[df$MethDiff <= -0.2 & df$log2FC >= 2] <- 'HypoMeth-UPExp'
df$status[df$MethDiff <= -0.2 & df$log2FC <= -2] <- 'HypoMeth-DownExp'
df$status[df$MethDiff >= 0.2 & df$log2FC >= 2] <- 'HyperMeth-UPExp'
df$status[df$MethDiff >= 0.2 & df$log2FC <= -2] <- 'HyperMeth-DownExp'
df$status[ov$subjectHits %in% ovPromoter$subjectHits & df$status != 'No-Sig'] <- 'SigPromoter'

message("Filtering out NA log2FC")
outlier <- boxplot(df$log2FC)
df <- df[!df$log2FC %in% outlier$out, ]


message("Plot Scatter Plot")
p <- ggplot(df, aes(x=MethDiff, y=log2FC, color=status)) + geom_point(size = 0.5) +
        scale_color_manual(values = c("#e07a5f", "#f2cc8f", "#81b29a", "#2d93ad", "#e0e0e0", '#3d405b')) +
        geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black") +
        xlab("Mean Methylation Beta Value Difference") + ylab("Gene Expression Log2 FoldChange") +
        theme_minimal(base_size = 20)

p_anno <- p + geom_text_repel(data = df[df$status == 'SigPromoter', ], aes(label = GeneName), size = 3.5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))

if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/CpGGeneScatterPlot.pdf", width = 14, height = 10)
print(p)
dev.off()


pdf("./Figure/CpGGeneScatterPlot_Labeled.pdf", width = 14, height = 10)
print(p_anno)
dev.off()

## Blow Part is for generating DMP-related Genes/Promoters

message("Prepare DMP related Gene File")
message("Prepare DMP (cutoff 15 meth.diff and q vlaue 0.05)")
myDiff15p <- getData(getMethylDiff(myDiff, difference=15, qvalue=0.05))
message("!Important, keey beta matrix filtered with DMP")
rownames(filteredBeta) <- rownames(myDiff)
beta15p <- filteredBeta[rownames(myDiff15p),]

message("Loading Annotation")
gff <- gffToGRanges("../3.DMP/Danio_rerio.GRCz11.98.gtf.gz")
grl <- as(split(gff, gff$type), "GRangesList")
grl[["promoter"]] <- promoters(grl[["gene"]])
source("../4.DMPVis/ExpandGeneRange.R")
expandedGenes <- ExpandGeneRange(grl[["gene"]])

message("Applying MatchingGene function to get DMP-related Genes/Promoters")
source("./MatchingGene.R")
DMP.ExpandedGene <- MatchingGene(myDiff15p, beta15p, expandedGenes)
DMP.Promoter <- MatchingGene(myDiff15p, beta15p, grl[["promoter"]])

AllProbe.ExpandedGene <- MatchingGene(getData(myDiff), filteredBeta, expandedGenes)
AllProbe.Promoter <- MatchingGene(getData(myDiff), filteredBeta, grl[["promoter"]])

if (!file.exists("./csv")) dir.create("./csv")
write.csv(DMP.ExpandedGene$Genes, "./csv/DMP_Related_ExpandedGenes.csv", quote=FALSE, row.names=FALSE)
write.csv(DMP.ExpandedGene$DMP, "./csv/DMP_Annotated_by_ExpandedGenes.csv", quote=FALSE, row.names=FALSE)

write.csv(DMP.Promoter$Genes, "./csv/DMP_Related_Promoters.csv", quote=FALSE, row.names=FALSE)
write.csv(DMP.Promoter$DMP, "./csv/DMP_Annotated_by_Promoter.csv", quote=FALSE, row.names=FALSE)

write.csv(AllProbe.ExpandedGene$Genes, "./csv/AllProbe_Related_ExpandedGenes.csv", quote=FALSE, row.names=FALSE)
write.csv(AllProbe.ExpandedGene$DMP, "./csv/AllProbe_Annotated_by_ExpandedGenes.csv", quote=FALSE, row.names=FALSE)

write.csv(AllProbe.Promoter$Genes, "./csv/AllProbe_Related_Promoters.csv", quote=FALSE, row.names=FALSE)
write.csv(AllProbe.Promoter$DMP, "./csv/AllProbe_Annotated_by_Promoter.csv", quote=FALSE, row.names=FALSE)

if (!file.exists("./Data")) dir.create("./Data")
save(DMP.ExpandedGene, file="./Data/DMP_ExpandedGene.rda")
save(DMP.Promoter, file="./Data/DMP_Promoter.rda")
save(AllProbe.ExpandedGene, file="./Data/AllProbe_ExpandedGene.rda")
save(AllProbe.Promoter, file="./Data/AllProbe_Promoter.rda")

message("Draw DMP-related ExpandedGene/Promoter Heatmap")
library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)

message("Plot at least 5 DMP enriched Promoter genes.")
dfGene <- DMP.Promoter$Genes
dfGene <- dfGene[dfGene$EnrichNumber >= 5, ]
geneMatrix <- as.matrix(dfGene$MeanBeta)
rownames(geneMatrix) <- dfGene$gene_name

annotation <- data.frame(pheno=rep(c("WildType", "Mutation"), each=3))
rownames(annotation) <- colnames(geneMatrix)
annoCol  <- c('#56B4E9', '#E69F00')
names(annoCol) <- c("WildType", "Mutation")
anno_colors <- list(pheno = annoCol)

colfunc <- colorRampPalette(c("#81b29a", "white", "#e07a5f"))

if (!file.exists("./Figure")) dir.create("./Figure")
graphics.off()
pdf("./Figure/AtLeast5DMPEnrichedPromoter.pdf", width = 6, height = 15)
pheatmap(geneMatrix, 
         color=colfunc(100), 
         annotation = annotation, 
         annotation_colors = anno_colors, 
         show_colnames=FALSE, 
         cluster_cols=FALSE)
dev.off()


dfGene <- DMP.ExpandedGene$Genes
dfGene <- dfGene[dfGene$EnrichNumber >= 5, ]
geneMatrix <- as.matrix(dfGene$MeanBeta)
rownames(geneMatrix) <- dfGene$gene_name

#labelName <- DMP.ExpandedGene$Genes$gene_name
#labelName[!labelName %in% c("vezf1a","pax2b", "epb41l4a", "sox5")] <- ""
#colfunc <- colorRampPalette(c("#81b29a", "white", "#e07a5f"))
#pheatmap(as.matrix(DMP.ExpandedGene$Genes$MeanBeta)[1:100, ], color=colfunc(100), main="DMP-related Genes")


message("A very nice script to plot Gene Methylation Heatmap Togather")

plot_list=list()
for (i in 1:5){
    trigger <- FALSE
    if(i == 5) trigger <- TRUE
    message("Start:", ((i-1) * 102 + 1), "   End:", (i * 102))
    x <- pheatmap(geneMatrix[((i-1) * 102 + 1) : (i * 102),],
                  color=colfunc(100),
                  annotation = annotation,
                  annotation_colors = anno_colors,
                  legend = trigger,
                  annotation_legend=trigger,
                  cluster_cols = FALSE, 
                  show_colnames=FALSE)
    plot_list[[i]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}
g <- grid.arrange(grobs=plot_list, layout_matrix = rbind(c(1,1,2,2,3,3,4,4,5,5,5),
                                                         c(1,1,2,2,3,3,4,4,5,5,5)))
if (!file.exists("./Figure")) dir.create("./Figure")
graphics.off()
ggsave("./Figure/AtLeast5DMPEnrichedExpandedGene.pdf", g, width=20, height=16)

