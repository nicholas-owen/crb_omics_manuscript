############################################################
# Aurthor : Tian
# Description: This script is written to show some global pattern between RNA-seq and Methylation
############################################################

library("methylKit")
library("GenomicRanges")
library("genomation")
library("ggplot2")
library("ggpubr")

message("Prepare TPM Matrix")
TPM <- read.csv("./External/TPM.tsv", sep="\t", as.is=T)
rownames(TPM) <- TPM$gene_id

message("Prepare DMP list")
load("../3.DMP/Data/myDiff.rda")
diff <- getData(myDiff)
# diff <- getData(getMethylDiff(myDiff, difference = 15, qvalue = 0.05))

message("Prepare Methylation Beta Value")
load("../3.DMP/Data/filteredBeta.rda")
beta <- filteredBeta

message("Loading Annotation")
gff <- gffToGRanges("../3.DMP/Danio_rerio.GRCz11.98.gtf.gz")
grl <- as(split(gff, gff$type), "GRangesList")

message("Mapping Location between Methylatio to Expression")
anno <- grl[['gene']]
message("Filter Gene Annotation with DEG list")
anno <- anno[as.numeric(seqnames(anno)) %in% c(1:25)]

message("Prepare ExtendedGene, Promoter, and Genebody three Annotation for mapping")
source("../4.DMPVis/ExpandGeneRange.R")
expandedGenes <- ExpandGeneRange(anno,  upstream = 2000)
Promoters <- promoters(anno)

message("Find overlap between ExtendedGene/Promoter/Genebody with Methylation")
ov <- suppressWarnings(as.data.frame(findOverlaps(expandedGenes, as(diff, "GRanges"))))
ovPromoter <- suppressWarnings(as.data.frame(findOverlaps(Promoters, as(diff, "GRanges"))))
ovGeneBody <- ov[!ov[,2] %in% ovPromoter[,2],]


methPheno <- c(0,0,0,1,1,1)
tpmPheno <- c(1,0,0,0,0,0,1,1,1,1,1)

message("Find Calculate Correlation for ExpendGene, Promoter, GeneBody")
source("./CorMeth.R")
source("./DrawRegression.R")
CorData <- CorMeth(beta, TPM, methPheno, tpmPheno, ov, expandedGenes)
ExpendGenePlot <- DrawRegression(CorData, "Gene + Promoter")


CorData <- CorMeth(beta, TPM, methPheno, tpmPheno, ovPromoter, expandedGenes)
PromoterPlot <- DrawRegression(CorData, "Promoter")

CorData <- CorMeth(beta, TPM, methPheno, tpmPheno, ovGeneBody, expandedGenes)
GeneBodyPlot <- DrawRegression(CorData, "GeneBody")

p <- ggarrange(ExpendGenePlot, PromoterPlot, GeneBodyPlot, ncol = 3, nrow = 1, align = "hv", widths = c(3, 3, 3), common.legend = TRUE)

if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/GlobalCorrelation.pdf", width = 12, height = 4)
print(p)
dev.off()

#### Below code is for plotting the TSS plot
message("Plot TSS Plot")

diff <- getData(myDiff)

# In below code I calcualted the DMP plot, but it's a bit pointless, as the meth.diff are always high.
# diff <- getData(getMethylDiff(myDiff, difference = 15, qvalue = 0.05))

message("Consider TSS as The Start Site of Transcripts")
TSS2000 <- promoters(grl[['transcript']], upstream=2000, downstream=2000)
ov <- suppressWarnings(as.data.frame(findOverlaps(TSS2000, as(diff, "GRanges"))))
ov$TSS <- (start(TSS2000[ov$queryHits,]) + end(TSS2000[ov$queryHits,]) - 1) / 2
ov$strand <- as.character(strand(TSS2000[ov$queryHits,]))
ov$CpGPosition <- diff$start[ov$subjectHits]
ov$dist <- ov$CpGPosition - ov$TSS
ov$meth.diff <- diff$meth.diff[ov$subjectHits]
ov$betaWT <- rowMeans(beta[ov$subjectHits, which(methPheno == 0)])
ov$betaMU <- rowMeans(beta[ov$subjectHits, which(methPheno == 1)])

message("Adjust CpG-to-TSS Distance based on Strand")
ov$dist[ov$strand == '-'] <- -1 * (ov$dist[ov$strand == '-'])

message("Group all CpG's into 1000 groups based on their distance to TSS")
groupInfo <- cut(ov$dist,1000, include.lowest=TRUE)

message("Plotting Meth Diff plot around TSS")
df <- data.frame(position=seq(-2000,2000,length.out=1000),
                 meth.diff=aggregate(ov$meth.diff, by=list(groupInfo), function(x) mean(x))$x,
                 betaWT=aggregate(ov$betaWT, by=list(groupInfo), function(x) mean(x))$x,
                 betaMU=aggregate(ov$betaMU, by=list(groupInfo), function(x) mean(x))$x)

p <- ggplot(df, aes(position, meth.diff)) +
  geom_point(color='grey') + 
  theme_minimal(base_size = 14) +
  labs(x = "Grouped of CpGs' distance to TSS") + 
  ylab("Averaged Methylation Diff (from methylKit)") +
  geom_smooth(span = 50)

if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/MethDiffAroundTSS.pdf", width = 12, height = 6)
print(p)
dev.off()

message("Plotting Beta Value plot around TSS")
df2 <- data.frame(position=c(df$position, df$position),
                  beta=c(df$betaWT, df$betaMU),
                  pheno=rep(c("WildType"," Mutation"), each=nrow(df)))

p <- ggplot(df2, aes(position, beta, colour=pheno)) +
  geom_point() + 
  theme_minimal(base_size = 14) +
  labs(x = "Grouped of CpGs' distance to TSS") + 
  ylab("Averaged Methylation Diff (from methylKit)") +
  geom_smooth(span = 50) +
  scale_color_manual(values=c('#E69F00', '#56B4E9'))

if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/BetaAroundTSS.pdf", width = 12, height = 6)
print(p)
dev.off()

message("Get Methylation Value Around CpG Island")

message("Load CGI information, forming Upstream-CGI, Downstream-CGI")
cgi <- read.csv("./External/cpgIslandExt.txt", header=F, sep="\t", as.is=T)[,c(2,3,4)]
colnames(cgi) <- c("seqname","start","end")
cgi$seqname <- substr(cgi$seqname,4,100)

cgi1 <- cgi2 <- cgi

cgi1$strand <- "-"
cgi1 <- makeGRangesFromDataFrame(cgi1)
cgi2$strand <- "+"
cgi2 <- makeGRangesFromDataFrame(cgi2)

cgi <- makeGRangesFromDataFrame(cgi)
cgiTail <- promoters(cgi1, upstream=2000, downstream=0)
cgiHead <- promoters(cgi2, upstream=2000, downstream=0)

message("Overlap 3-stage CGI to myDiff (all CpGs)")
CGI <- suppressWarnings(as.data.frame(findOverlaps(cgi, as(diff, "GRanges"))))
Head <- suppressWarnings(as.data.frame(findOverlaps(cgiHead, as(diff, "GRanges"))))
Tail <- suppressWarnings(as.data.frame(findOverlaps(cgiTail, as(diff, "GRanges"))))

message("Read function to calculate Average value for 3-stage CGI")
GetSmoothValue <- function(cgi,Body, mydiff)
{
    Body$CGI <- (start(cgi[Body[,1],]) + end(cgi[Body[,1],]))/2
    Body$CGLength <- abs(start(cgi[Body[,1],])  - (end(cgi[Body[,1],])))
    Body$CpGPosition <- mydiff$start[Body[,2]]
    Body$dist <- (Body$CpGPosition - Body$CGI)/(Body$CGLength/2)
    Body$meth.diff <- mydiff$meth.diff[Body[,2]]

    beta_WT <- rowMeans(beta[,1:3])
    beta_Mutation <- rowMeans(beta[,4:6])

    groupInfo_CGI <- cut(Body$dist,400,include.lowest=TRUE)

    MeanMeth_CGI <- aggregate(Body$meth.diff, by=list(groupInfo_CGI),function(x) mean(x))
    MeanMeth_WT <- aggregate(beta_WT[Body[,2]],by=list(groupInfo_CGI),function(x) mean(x))
    MeanMeth_Mutation <- aggregate(beta_Mutation[Body[,2]],by=list(groupInfo_CGI),function(x) mean(x))

    return(list(MethDiff=MeanMeth_CGI$x, BetaWT=MeanMeth_WT$x, BetaMutation=MeanMeth_Mutation$x))
}

message("Calculate 3-stage Value")
HeadLine <- GetSmoothValue(cgiHead, Head, diff)
BodyLine <- GetSmoothValue(cgi, CGI, diff)
TailLine <- GetSmoothValue(cgiTail, Tail, diff)

df3 <- data.frame(position=1:1200,
                  MethDiff=c(HeadLine$MethDiff, BodyLine$MethDiff, TailLine$MethDiff),
                  BetaWT=c(HeadLine$BetaWT, BodyLine$BetaWT, TailLine$BetaWT),
                  BetaMU=c(HeadLine$BetaMutation, BodyLine$BetaMutation, TailLine$BetaMutation))

p <- ggplot(df3, aes(position, MethDiff)) +
  geom_point(color='grey') +
  theme_minimal(base_size = 14) +
  labs(x = "Grouped of CpGs' distance to CpG Islands") +
  ylab("Averaged Methylation Diff (from methylKit)") +
  geom_smooth() + 
  scale_x_continuous(breaks = c(1, 400, 600, 800, 1150), labels = c("UpStream 2000", "|", "CpG Islands", "|", "DownStream 2000")) + 
  geom_vline(xintercept = c(400, 800), linetype = "dashed", color = "grey")

if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/MethDiffAroundCpGIsland.pdf", width = 12, height = 6)
print(p)
dev.off()


message("Plotting Beta value across CpG Islands")
df4 <- data.frame(position=rep(1:1200, 2),
                  beta=c(df3$BetaWT, df3$BetaMU),
                  pheno=rep(c("WildType"," Mutation"), each=nrow(df3)))

p <- ggplot(df4, aes(position, beta, colour=pheno)) +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(x = "Grouped of CpGs' distance to CpG Islands") +
  ylab("Averaged Beta") +
  geom_smooth() +
  scale_color_manual(values=c('#E69F00', '#56B4E9')) +
  scale_x_continuous(breaks = c(1, 400, 600, 800, 1200), labels = c("UpStream 2000", "|", "CpG Islands", "|", "DownStream 2000")) +
  geom_vline(xintercept = c(400, 800), linetype = "dashed", color = "grey")


if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/BetaAroundCpGIslands.pdf", width = 12, height = 6)
print(p)
dev.off()

