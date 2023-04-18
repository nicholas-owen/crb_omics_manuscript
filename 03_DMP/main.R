#############################################################
# Description: A Script for calculation of Differential 
# Methylation Probes and Differential Methylation Regions
#############################################################

library("methylKit")

message("Removed all CpGs have 1 or 0 value across all Sampels")
load("../02_Methylkit/Data/meth.rda")
load("../02_Methylkit/Data/beta.rda")

index <- which(!rowSums(beta) %in% c(0,ncol(beta)))
meth.dmp <- meth[index]
myDiff <- calculateDiffMeth(meth.dmp, overdispersion="MN", mc.cores=40)

message("Draw P value Distribution")
if (!file.exists("./Figure")) dir.create("./Figure")
pdf("./Figure/PvalueDistribution.pdf", width=8, height=8)
hist(myDiff$pvalue,border=0, main="P value Distribution", xlab=paste(nrow(myDiff),'CpGs'))
dev.off()

message("Filtering DMPs across all CpGs")
myDiff15p <- getMethylDiff(myDiff, difference=15, qvalue=0.05)

message("Save Filtered Beta and Meth")
if (!file.exists("./Data")) dir.create("./Data")
filteredBeta <- beta[index,]
filteredMeth <- meth[index]
save(filteredBeta, file="./Data/filteredBeta.rda")
save(filteredMeth, file="./Data/filteredMeth.rda")
save(myDiff,file="./Data/myDiff.rda")

message("Annotate DMPs")
library("GenomicRanges")
library("genomation")
system("wget ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz")
gff <- gffToGRanges("./Danio_rerio.GRCz11.98.gtf.gz")
grl <- as(split(gff, gff$type), "GRangesList")

tmpDiff <- getData(myDiff)
ov <- findOverlaps(as(tmpDiff, "GRanges"), grl[["gene"]])
Anno <- as.data.frame((grl[['gene']][as.data.frame(ov)$subjectHits ,c("gene_id", "gene_name")]))

tmpDiff$ovIndex <- 1:nrow(tmpDiff)
Anno$ovIndex <- as.data.frame(ov)$queryHits
AnnotatedDMP <- merge(tmpDiff, Anno, by="ovIndex", all = T)
colnames(AnnotatedDMP)[c(2:5, 9:12)] <- c("cpg_chr", "cpg_start", "cpg_end", "cpg_strand", "gene_chr", "gene_start", "gene_end", "gene_strand")

message("Save All CpGs and Annotation Into CSV")
if (!file.exists("./csv")) dir.create("./csv")
write.csv(AnnotatedDMP, "./csv/AnnotatedDMP.csv", quote=FALSE, row.names=FALSE)

