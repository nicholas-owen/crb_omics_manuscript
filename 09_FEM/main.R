############################################################
# Description: An implementation of RRBS version FEM
############################################################

library("methylKit")
library("genomation")

message("Prepare RNA-seq Results")
DEG <- read.csv("../06_IntergrateRNAseq/External/deseq2-HF-M56-annotated.diffexpr-results.csv", as.is=TRUE)
rownames(DEG) <- DEG$gene_id

message("Prepare PPI Network")
load("./Data/ppi.rda")

message("Prepare DMP-related Genes")
load("../06_IntergrateRNAseq/Data/AllProbe_ExpandedGene.rda")
load("../06_IntergrateRNAseq/Data/AllProbe_Promoter.rda")

DMPGeneList <- list(AllProbe.ExpandedGene, AllProbe.Promoter)
pheno.v <- c("WT", "WT", "WT", "crb2a", "crb2a" ,"crb2a")
names(DMPGeneList) <- c("AllProbeExpandedGene", "AllProbePromoter")

source("./RunRRBSFEM.R")

if (!file.exists("./Data")) dir.create("./Data")

break

for(i in names(DMPGeneList))
{
    message("Calculating ", i)
    FEMResult <- RunRRBSFEM(DMPGeneList[[i]], pheno.v)
    save(FEMResult, file=paste0("./Data/", i , '_FEMResult.rda'))
}

for(i in names(DMPGeneList))
{
    library("marray")
    library("corrplot")
    source("./myFemModShow.R")
    
    if (!file.exists("./Figure")) dir.create("./Figure")
    if (!file.exists(paste0("./Figure/",i))) dir.create(paste0("./Figure/",i))
    load(paste0("./Data/", i, "_FEMResult.rda"))
    setwd(paste0("./Figure/",i,"/"))
    for(i in names(FEMResult$topmod)) {
        FemModShow(FEMResult$topmod[[i]], name=i, FEMResult)
        fwrite(FEMResult$topmod[[i]], paste0("FEM_top_result_promoter", i, "_table.csv"))
    }
    setwd("../../")
}


# Gene specific analysis
# - cdk6
gnel<-FemModShow(FEMResult$topmod[['cdk6']], name="cdk6", FEMResult)
ig<-graph_from_graphnel(gnel, name = TRUE, weight = TRUE, unlist.attrs =  TRUE)
BiocManager::install("NetPathMiner")
library("NetPathMiner")
plotNetwork(ig, vertex.color="compartment.name")
plotCytoscapeGML(ig, file="example.gml", layout = layout.c,
                 vertex.size = 5, vertex.color = v.color)
# -cdk6
library(RCy3)
createNetworkFromGraph(gnel,"myGraph4")
# merge network table at:./Figure/AllProbeExpandedGene/csv/FEM-cdk6-cyto_table.csv with exp table
df.network <- read.csv("./Figure/AllProbeExpandedGene/csv/FEM-cdk6-cyto_table.csv")
df.exp<- read.csv("../06_IntergrateRNAseq/External/deseq2-HF-M56-annotated.diffexpr-results.csv", as.is=TRUE)
df.exp<-df.exp%>%
    dplyr::select(c("gene_id", "gene_name", "log2FoldChange","padj"))
colnames(df.exp)[1]<-"id"
df.netexp<-dplyr::left_join(df.network, df.exp)
write.csv(df.netexp, file="./Figure/AllProbeExpandedGene/csv/FEM-cdk6-cyto_table_annot.csv",quote = FALSE, row.names = FALSE)

# - bmpr1aa
gnel<-FemModShow(FEMResult$topmod[['bmpr1aa']], name="bmpr1aa", FEMResult)
createNetworkFromGraph(gnel,"bmpr1aa")
df.network <- read.csv("./Figure/AllProbeExpandedGene/csv/FEM-bmpr1aa-cyto_table.csv")
df.network<-df.network %>%
    dplyr::select(-gene_name)
df.exp<- read.csv("../06_IntergrateRNAseq/External/deseq2-HF-M56-annotated.diffexpr-results.csv", as.is=TRUE)
df.exp<-df.exp%>%
    dplyr::select(c("gene_id", "gene_name", "log2FoldChange","padj"))
colnames(df.exp)[1]<-"id"
df.netexp<-dplyr::left_join(df.network, df.exp, by="id")
write.csv(df.netexp, file="./Figure/AllProbeExpandedGene/csv/FEM-bmpr1aa-cyto_table_annot.csv",quote = FALSE, row.names = FALSE)

# - gstk2
gnel<-FemModShow(FEMResult$topmod[['gstk2']], name="gstk2", FEMResult)
createNetworkFromGraph(gnel,"gstk2")
df.network <- read.csv("./Figure/AllProbeExpandedGene/csv/FEM-gstk2-cyto_table.csv")
df.network<-df.network %>%
    dplyr::select(-gene_name)
df.exp<- read.csv("../6.IntergrateRNAseq/External/deseq2-HF-M56-annotated.diffexpr-results.csv", as.is=TRUE)
df.exp<-df.exp%>%
    dplyr::select(c("gene_id", "gene_name", "log2FoldChange","padj"))
colnames(df.exp)[1]<-"id"
df.netexp<-dplyr::left_join(df.network, df.exp, by="id")
write.csv(df.netexp, file="./Figure/AllProbeExpandedGene/csv//FEM-gstk2-cyto_table_annot.csv",quote = FALSE, row.names = FALSE)



