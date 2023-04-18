# 03 Differentially Methylated Positions (DMP) Analysis

DMPs were called using the `methylkit` software. Note that we discarded CpGs with all 0 or all 1 methylation status. _P_ value is adjusted for multiple testing with _q_ value.

The Ensembl GRCz11 annotation was used from [ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz](ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz) and with `genomation` R package.

## CSV folder:

 - `AnnotatedDMP.csv`: fully annotated DMP information, contains 717343 lines. It included all CpGs (DMP or not), with their corresponding gene (if mapped).

## Data folder:

- `filteredBeta.rda`: beta matrix after filtering out uninformative CpGs.

- `filteredMeth.rda`: meth matrix after filtering out uninformative CpGs.

- `myDiff.rda`: Differential Calling results, contains all CpGs, their significance _p_ value, meth.diff, etc.

## Figure folder:

- `PvalueDistribution.pdf`: The _p_ value distribution from all CpGs.
