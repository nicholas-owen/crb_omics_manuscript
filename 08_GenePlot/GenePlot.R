
library("gridExtra")
library("ggplot2")
library("ggpubr")

GenePlot <- function(Diff, beta, TPM, DMG, anno, GeneID) {
    x <- GeneID
    GeneName <- DMG[x, "gene_name"]
    message("Printing ", x, " | ", GeneName, " ...")
    
    drawData <- data.frame(tpm = as.numeric(TPM[x, 1:11]), pheno = tpmPD)
    
    fill <- "white"
    line <- "#1F3552"
    
    tpmPlot <- ggplot(drawData, aes(x = pheno, y = tpm, col = pheno)) + geom_boxplot(fill = fill, 
        alpha = 0.7, outlier.shape = 20) + geom_jitter(width = 0.2, size = 2) + labs(x = "", 
        y = "Expression TPM Value") + ggtitle(GeneName) + theme(plot.title = element_text(hjust = 0.5)) + 
        theme_bw() + annotate("text", x = 1.5, y = mean(range(drawData$tpm)), label = paste("P-value = ", 
        WriteOutPval(DEG[x, "padj"]), sep = "")) + scale_color_manual(values = c("#E69F00", 
        "#56B4E9"))
    
    ol <- suppressWarnings(as.data.frame(findOverlaps(anno[which(anno$gene_id == 
        x)], as(Diff, "GRanges"))))
    
    drawData <- data.frame(beta = as.numeric(beta[ol[, 2], ]), position = 1:nrow(ol), 
        pheno = c(rep("WildType", nrow(ol) * 3), rep("Mutation", nrow(ol) * 3)))
    
    scatterPlot <- ggplot(drawData, aes(x = position, y = beta, color = pheno)) + 
        geom_point(size = 1.8, alpha = 0.8) + labs(x = "", y = "methylation beta value") + 
        scale_x_continuous(breaks = c(1:nrow(ol)), labels = myDiff[ol[, 2]]$start) + 
        geom_smooth(method = loess) + ggtitle(GeneName) + theme(plot.title = element_text(hjust = 0.5)) + 
        theme_bw() + annotate("text", x = nrow(ol)/2, y = mean(range(drawData$beta)), 
        label = paste("Meth.Diff = ", round(DMG[x, "MeanMethDiff"], 2), sep = "")) + 
        scale_color_manual(values = c("#E69F00", "#56B4E9"))
    
    yplot <- ggdensity(drawData, "beta", fill = "pheno", palette = c("#E69F00", "#56B4E9")) + 
        rotate()
    
    print(ggarrange(tpmPlot, scatterPlot, yplot, ncol = 3, nrow = 1, align = "hv", 
        widths = c(3, 8, 2), common.legend = TRUE))
}
