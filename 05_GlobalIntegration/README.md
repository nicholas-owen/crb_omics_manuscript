# 05 Global Integration

This is the first work to integrate CpG information from RRBS with RNA-seq and genome. It will calculate global correlation across all CpG sites and related gene expression `TPM` data and plot CpG methylation different status around TSS and CpG Islands. Finally it will plot the CpG methylation beta value around TSS and CpG Islands.

- `CorMeth.R`: The function to calculate correlation plot value from `TPM` and `beta` value.

- `DrawRegression.R`: Plot the two line correlation plot from results of `CorMeth.R`

## External folder:

Files generated externally.

- `cpgIslandExt.txt:` CpG Island information downloaded from `USCS`.

- `TPM.tsv`: `TPM` matrix provided from the bulk RNA-seq data analysis.

## Figure Folder:

- `BetaAroundCpGIslands.pdf`: `Beta` value group CpG island (upstream and downstream 2000).

- `GlobalCorrelation.pdf`: Correlation plot from all CpG's methylation and all gene's expression `TPM`.

- `MethDiffAroundTSS.pdf`: Methylation difference around TSS. TSS is calcualted from `genomation` package `transcript` list.

- `BetaAroundTSS.pdf`: `Beta` value around TSS.

- `MethDiffAroundCpGIsland.pdf`: Methylation difference around CpG island (upstream and downstream 2000).
