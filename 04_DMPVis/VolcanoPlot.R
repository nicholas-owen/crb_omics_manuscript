############################################################
# Description: A Script for Volcano plot
############################################################

VolcanoPlot <- function(diff, gff) {
    source("./ExpandGeneRange.R")
    
    library("ggplot2")
    library("ggrepel")
    
    grl <- as(split(gff, gff$type), "GRangesList")
    expandedGenes <- ExpandGeneRange(grl[["gene"]])
    ov <- suppressWarnings(as.data.frame(findOverlaps(expandedGenes, makeGRangesFromDataFrame(diff))))
    
    df <- diff
    df$stats <- "No-Sig"
    df$stats[df$qvalue <= 0.05 & df$meth.diff >= 15] <- "Hyper-Meth"
    df$stats[df$qvalue <= 0.05 & df$meth.diff <= -15] <- "Hypo-Meth"
    df$gene <- ""
    df$gene[ov[, 2]] <- as.character(as.data.frame(expandedGenes)[ov[, 1], "gene_name"])
    
    p <- ggplot(df, aes(x = meth.diff, y = -log10(qvalue), col = stats)) + geom_point(size = 0.5) + 
        scale_color_manual(values = c("#e07a5f", "#81b29a", "#999999")) + geom_hline(yintercept = c(-log10(0.05), 
        100), linetype = "dashed", color = "black") + geom_vline(xintercept = c(-15, 
        15), linetype = "dashed", color = "black") + theme_minimal(base_size = 14) + 
        geom_text_repel(data = df[-log10(df$qvalue) >= 100, ], aes(label = gene), 
            size = 3.5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, 
                "lines"))
    
    if (!file.exists("./Figure")) 
        dir.create("./Figure")
    pdf("./Figure/DMPVolcanoPlot.pdf", width = 10, height = 10)
    print(p)
    dev.off()
}
