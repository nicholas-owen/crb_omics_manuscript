############################################################ 
# Description: Differential Methylated and Differential 
# Expression Genes. Plot Venn.
############################################################.

library("methylKit")
library("GenomicRanges")
library("genomation")
library("ggplot2")
library("ggrepel")

message("Prepare RNA-seq Results")
Expression <- read.csv("../06_IntergrateRNAseq/External/deseq2-HF-M56-annotated.diffexpr-results.csv", 
    as.is = TRUE)
rownames(Expression) <- Expression$gene_id

message("Load Methylation Related Genes")
load("../06_IntergrateRNAseq/Data/DMP_ExpandedGene.rda")
load("../06_IntergrateRNAseq/Data/DMP_Promoter.rda")

message("Prepare TPM Matrix")
TPM <- read.csv("../05_GlobalIntegration/External/TPM.tsv", sep = "\t", as.is = T)
rownames(TPM) <- TPM$gene_id

message("Get Differential Expression Genes, with padj: 0.05, |log2FC| >= 2, baseMean >=10")
DE <- Expression[Expression$padj <= 0.05 & abs(Expression$log2FoldChange) >= 2 & 
    Expression$baseMean >= 10, ]
message("Get Differential Methylated Genes, with at least 5 CpGs Enriched")
DM <- DMP.ExpandedGene$Genes  # This list could be changed into any results from MatchGene function in ../6.IntergrateRNAseq


message("Finding DEDM Genes.")
DM <- DM
DEDM <- intersect(rownames(DE), rownames(DM))
DEDM <- cbind(DE[DEDM, ], DM[DEDM, ])

message("Save DEDM Genes and Annotation Into CSV")
if (!file.exists("./csv")) dir.create("./csv")
write.csv(DEDM, "./csv/DEDM.csv", quote = FALSE, row.names = FALSE)

message("Finding Strict DEDM Genes")
StrictDM <- DM[DM$EnrichNumber >= 5, ]
StrictDEDM <- intersect(rownames(DE), rownames(StrictDM))
StrictDEDM <- cbind(DE[StrictDEDM, ], StrictDM[StrictDEDM, ])

message("Save Stricr DEDM CpGs and Annotation Into CSV")
if (!file.exists("./csv")) dir.create("./csv")
write.csv(StrictDEDM, "./csv/StrictDEDM.csv", quote = FALSE, row.names = FALSE)

message("Plot DEDM and Strict DEDM Genes")
source("./PlotDEDMHeatmap.R")

g <- PlotDEDMHeatMap(StrictDEDM, TPM)
g_label <- PlotDEDMHeatMap(StrictDEDM, TPM, c("vsnl1b", "epb41l4a", "pcdh1gb9")) ### Please change gene name into whatever you like

if (!file.exists("./Figure")) dir.create("./Figure")
graphics.off()
ggsave("./Figure/StrictDEDMHeatmap_full_label.pdf", g, width=8.5, height=8.5)
ggsave("./Figure/StrictDEDMHeatmap_custom_label.pdf", g, width=8.5, height=8.5)


g_dedm_label <- PlotDEDMHeatMap(DEDM, TPM, c("vezf1a", "pax2b", "epb41l4a", "MAP3K13")) ### Please change gene name into whatever you like
if (!file.exists("./Figure")) dir.create("./Figure")
graphics.off()
ggsave("./Figure/DEDMHeatmap_custom_label.pdf", g_dedm_label, width=8.5, height=8.5)

message("Plot Venn Plot")

Hyper <- DM$gene_id[DM$MeanMethDiff > 0]
Hypo <- DM$gene_id[DM$MeanMethDiff < 0]
UpRegulate <- DE[which(DE$log2FoldChange >= 2 & DE$padj <= 0.05), "gene_id"]
DownRegulate <- DE[which(DE$log2FoldChange <= -2 & DE$padj <= 0.05), "gene_id"]

library(VennDiagram)
HyperPlotVenn <- venn.diagram(
  x = list(Hyper, UpRegulate, DownRegulate),
  category.names = c("HyperMeth" , "UpExp", "DownExp"),
  resolution = 300,
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  cex = 2,
  fontface = "bold",
  fill = c("#e07a5f", "#f40000", "#00ff00"),
  cat.cex = 1,
  cat.fontface = "bold",
  filename=NULL)

library(grDevices)
if (!file.exists("./Figure")) dir.create("./Figure")
graphics.off()
pdf(file="./Figure/HyperVennPlot.pdf")
    grid.draw(HyperPlotVenn)
dev.off()


HypoPlotVenn <- venn.diagram(
  x = list(Hypo, UpRegulate, DownRegulate),
  category.names = c("HypoMeth" , "UpExp", "DownExp"),
  resolution = 300,
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  cex = 2,
  fontface = "bold",
  fill = c("#81b29a", "#f40000", "#00ff00"),
  cat.cex = 1,
  cat.fontface = "bold",
  filename=NULL)

library(grDevices)
if (!file.exists("./Figure")) dir.create("./Figure")
graphics.off()
pdf(file="./Figure/HypoVennPlot.pdf")
    grid.draw(HypoPlotVenn)
dev.off()

