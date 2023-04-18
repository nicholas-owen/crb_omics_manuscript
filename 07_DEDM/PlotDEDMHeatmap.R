
PlotDEDMHeatMap <- function(DEDM, TPM, labelName=NULL) {

    label <- DEDM$gene_name
    if(!is.null(labelName))
    {
        label[!label %in% labelName] <- ""
    }
    
    dfTPM <- as.matrix(TPM[rownames(DEDM), c(c(2:6), c(1, 7:11))])
    dfMeth <- as.matrix(DEDM$MeanBeta)
    rownames(dfTPM) <- rownames(dfMeth) <- DEDM$gene_name
    
    library("grid")
    library("gridExtra")
    library("pheatmap")
    library("ggplot2")
    
    message("Prepare for Methylation Heatmap Plotting")
    annotation <- data.frame(pheno = rep(c("WildType", "Mutation"), each = 3))
    rownames(annotation) <- colnames(dfMeth)
    annoCol <- c("#56B4E9", "#E69F00")
    names(annoCol) <- c("WildType", "Mutation")
    anno_colors <- list(pheno = annoCol)
    
    colfuncMeth <- colorRampPalette(c("#81b29a", "white", "#e07a5f"))
    colfuncTPM <- colorRampPalette(c(c("#00ff00", "black", "#f40000")))
    
    message("Plotting Pheatmap for Methylation Beta Matrix")
    tmp <- hclust(dist(dfMeth))
    pMeth <- pheatmap(dfMeth[tmp$order, ], color = colfuncMeth(100), annotation = annotation, 
        annotation_colors = anno_colors, main = "    Averaged Beta Value", annotation_legend = FALSE, 
        show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_row=label[tmp$order])
    
    
    box <- boxplot(as.vector(dfTPM))
    outer <- which(dfTPM %in% box$out)
    dfTPM[outer] <- max(dfTPM[-outer])
    
    annotation <- data.frame(pheno = c(rep("WildType", 5), rep("Mutation", 6)))
    rownames(annotation) <- colnames(dfTPM)
    annoCol <- c("#56B4E9", "#E69F00")
    names(annoCol) <- c("WildType", "Mutation")
    anno_colors <- list(pheno = annoCol)
    
    message("Plotting Pheatmap for Methylation TPM Matrix")
    pTPM <- pheatmap(dfTPM[tmp$order, ], color = colfuncTPM(100), annotation = annotation, 
        annotation_colors = anno_colors, main = "Gene Expression TPM Value", show_colnames = FALSE, 
        cluster_cols = FALSE, cluster_rows = FALSE, labels_row=label[tmp$order])
    
    g <- grid.arrange(grobs = list(pMeth[[4]], pTPM[[4]]), layout_matrix = rbind(c(1, 
        1, 2, 2, 2), c(1, 1, 2, 2, 2)))
    return(g)
}
