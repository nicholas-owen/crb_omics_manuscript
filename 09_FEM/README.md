# 09 FEM Analysis

This is script for doing `FEM` calculation on `RRBS` and `RNAseq` data. For algorithm detail, please read:

__A systems-level integrative framework for genome-wide DNA methylation and gene expression data identifies differential gene expression modules under epigenetic control, Jiao __et al__ 2014__ [https://pubmed.ncbi.nlm.nih.gov/24794928/](https://pubmed.ncbi.nlm.nih.gov/24794928/)


- `GenStatM.R:` The function to generate statistic results for methylation gene.

- `PreparePPI.R`: A script to prepare PPI network for all species.

- `RunRRBSFEM.R`: A function to run RRBS, take DMP-related gene list as input.

- `myDoFEMbi.R`: Function odified from FEM, to calculate FEM module.

- `myFemModShow.R`: Function modified from FEM, to plot modules.

## Figure folder:

- `Figure/AllProbeExpandedGene`: All plot for modules generate by all CpG mapped expanded gene list.

- `Figure/AllProbePromoter`: All plot for modules generate by all CpG mapped promoter list.

## Cytoscape folder

- `FEM-bmpr1aa.cys`: FEM hub network for _bmpr1b_ hub gene network

- `FEM-cdk6.cys`: FEM hub network for _cdk6_ hub gene network

Files are for [Cytoscape 3.8 +](http://www.cytoscape.org/) 