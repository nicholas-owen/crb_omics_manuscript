# CRB Omics Manuscript Repo
Code for the _crb2a_ multi omics manuscript

[![DOI](https://zenodo.org/badge/380179339.svg)](https://zenodo.org/badge/latestdoi/380179339)

---
This is the analysis pipeline for the DNA methylation analysis and integration of RNA-seq and methylation (RRBS) data for the manuscript :

### Loss of the crumbs cell polarity complex disrupts epigenetic transcriptional control and cell cycle progression in the developing retina

Owen, N., Toms, M., Tian, Y., Toualbi, L., Richardson, R., Young, R., Tracey-White, D., Dhami, P.,  Beck, S., Moosajee, M.

__Journal of Pathology 2023__ Apr;259(4):441-454. doi: [10.1002/path.6056](https://doi.org/10.1002/path.6056)

---

## RRBS data analysis pipeline


 The code is organised for reviewers and reseachers to validate our analysis result. The RRBS data was generated from control (wt) and _crb2a_<sup>-/-</sup> zebrafish at 56 hpf. The raw NuGen generated `bcl` file are taken as input. The process of data analysis has been outlined in the following numbered section, indicating the sequence of analysis. Input for subsequent sections often requires results from the formal analysis.

In each folder, there is a R script called `main.R`, by running it in R session, in theory all results presented in this repo could be reproduced. Users need to install corresponding R packages, install corresponding R version (4.0.2+).

All analysis within are using Ensembl genome as genomic annotation ([version 98](ftp://ftp.ensembl.org/pub/release-98/)), except for CpG islands information used in ./05_GlobalIntegration, which is downloaded from USCS.




## Data

The raw data for both RNA-seq and RRBS has been depositied within the NCBI Gene Expression Omnibus (GEO) at [GSE178709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178709) and [GSE178842](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178842)



