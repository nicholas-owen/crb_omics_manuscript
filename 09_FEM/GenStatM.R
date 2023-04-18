# GenStatM Function modified from FEM package, to calcualte RRBS.
library("limma")

GenStatM <- function(data.m, pheno.v){

    sampletype.f <- as.factor(pheno.v)

    design.sample <- model.matrix(~0 + sampletype.f)
    colnames(design.sample) <- levels(sampletype.f)
    sampletypes.v <- levels(sampletype.f)
    lmf.o <- lmFit(data.m, design.sample)
    ntypes <- length(levels(sampletype.f))
    ncomp <- 0.5 * ntypes * (ntypes - 1)
    cont.m <- matrix(0, nrow = ncol(design.sample), ncol = ncomp)
    tmp.v <- vector()

        c <- 1
    for (i1 in 1:(ntypes - 1)) {
        for (i2 in (i1 + 1):ntypes) {
            cont.m[i1, c] <- -1
            cont.m[i2, c] <- 1
            tmp.v[c] <- paste(sampletypes.v[i2], "--", sampletypes.v[i1],
                sep = "")
            c <- c + 1
        }
    }
    rownames(cont.m) <- sampletypes.v
    colnames(cont.m) <- tmp.v
    lmf2.o <- contrasts.fit(lmf.o, cont.m)
    bay.o <- eBayes(lmf2.o)
    top.lm <- list()

        for (c in 1:ncol(cont.m)) {
        top.lm[[c]] <- topTable(bay.o, coef = c, adjust.method = "fdr",
            number = nrow(data.m))
    }
        return(list(top = top.lm, cont = cont.m))

}
