############################################################
# Description: A Script for Visualise DMP
############################################################

library("methylKit")
library("GenomicRanges")
library("genomation")
library("ggplot2")

GeneFeatureBP <- function(diff, gff) {
    grl <- as(split(gff, gff$type), "GRangesList")
    grl[["promoter"]] <- promoters(grl[["gene"]])
    gff <- unlist(grl)  # Add Promoter into gff file
    
    hyperDMP <- diff[diff$meth.diff >= 0, ]
    hypoDMP <- diff[diff$meth.diff < 0, ]
    
    hyperNumber <- (suppressWarnings(annotateWithFeatures(as(hyperDMP, "GRanges"), 
        grl)))@num.annotation
    hypoNumber <- (suppressWarnings(annotateWithFeatures(as(hypoDMP, "GRanges"), 
        grl)))@num.annotation
    
    df <- data.frame(DMPstatus = rep(c("hyper", "hypo"), each = length(hyperNumber)), 
        feature = rep(names(hyperNumber), 2), number = c(hyperNumber, hypoNumber))
    
    if (!file.exists("./Figure")) 
        dir.create("./Figure")
    p <- ggplot(data = df, aes(x = feature, y = number, fill = DMPstatus)) + geom_bar(stat = "identity") + 
        theme_minimal(base_size = 14)
    # theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))
    
    p <- p + scale_fill_manual(values = c("#e07a5f", "#81b29a"))
    pdf("./Figure/GeneFeatureBarplot.pdf", width = 16, height = 8)
    print(p)
    dev.off()
}

