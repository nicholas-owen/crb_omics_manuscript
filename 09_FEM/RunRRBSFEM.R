
RunRRBSFEM <- function(DMPGenes, pheno.v)
{
    data.m <- DMPGenes$Genes$MeanBeta
    rownames(data.m) <- rownames(DMPGenes$Genes)
    source("./GenStatM.R")

    SM <- GenStatM(data.m, pheno.v)

    message("Prepare data for FEM function.")
    statM <- SM$top[[1]][,c('t', 'P.Value')]
    statR <- DEG[,c('stat','pvalue')]

    commonGene <- intersect(rownames(statM), rownames(ppi))
    commonGene <- commonGene[which(!is.na(statM[commonGene,1]) & !is.na(statR[commonGene,1]))]

    library('FEM')
    intFEM.o <- list(statM=statM[commonGene, ], statR=statR[commonGene, ],adj=ppi[commonGene,commonGene])
    source("./myDoFEMbi.R")
    FEMResult <- DoFEMbi(intFEM.o, nseeds=100,gamma=0.5,nMC=1000,sizeR.v=c(1,100), minsizeOUT=10,writeOUT=TRUE,nameSTUDY="Zebrafish-    Crb2a",ew.v=NULL);

    return(FEMResult)
}

