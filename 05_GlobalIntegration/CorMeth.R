
CorMeth <- function(beta, TPM, methPheno, tpmPheno, ov, anno)
{
    Mean.beta <- aggregate(beta[ov$subjectHits, ], by=list(ov$queryHits), function(x) mean(x))
    Mean.beta <- data.frame(Mean.beta, 
                            WildType_Meth=rowMeans(Mean.beta[, (which(methPheno == 0) + 1)]),
                            Mutation_Meth=rowMeans(Mean.beta[, (which(methPheno == 1) + 1)])) 

    GeneMeth <- data.frame(as.data.frame(anno[Mean.beta[,1], c("gene_id", "gene_name")]), Mean.beta)
    rownames(GeneMeth) <- GeneMeth$gene_id

    Mean.tpm <- data.frame(TPM, 
                           WildType_TPM=rowMeans(TPM[, which(tpmPheno == 0)]),
                           Mutation_TPM=rowMeans(TPM[, which(tpmPheno == 1)]))

    CommonGene <- intersect(rownames(Mean.tpm), rownames(GeneMeth))

    tmp <- data.frame(Mean.tpm[CommonGene,c("Mutation_TPM", "WildType_TPM")], GeneMeth[CommonGene,c("Mutation_Meth", "WildType_Meth")])

    q <- boxplot(tmp$Mutation_TPM)
    outlier <- which(tmp$Mutation_TPM %in% q$out)
    stmp <- tmp[-outlier,]
    groupInfo <- cut(stmp$Mutation_TPM, 30, include.lowest=TRUE)
    MeanMutation_Meth <- aggregate(stmp[,"Mutation_Meth"], by=list(groupInfo), function(x) mean(x))

    q <- boxplot(tmp$WildType_TPM)
    outlier <- which(tmp$WildType_TPM %in% q$out)
    stmp <- tmp[-outlier,]
    groupInfo <- cut(stmp$WildType_TPM, 30, include.lowest=TRUE)
    MeanWildType_Meth <- aggregate(stmp[,"WildType_Meth"], by=list(groupInfo), function(x) mean(x))

    return(list(MutationMeth=MeanMutation_Meth, WildTypeMeth=MeanWildType_Meth))
}

