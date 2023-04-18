# 06 Intergratation with RNA-seq data

Integrating analysis between differential methylation and differential expression. The key approach was to map methylation data to gene features, and generated gene's methylation value, for use in following analysis.

## CSV folder:

- `AllProbe_Annotated_by_ExpandedGenes.csv`: All probe annotated expanded gene's information.

- `DMP_Annotated_by_ExpandedGenes.csv`: DMP annotated expanded gene's information.

- `AllProbe_Annotated_by_Promoter.csv`: All probe annotated promoter's information.

- `DMP_Annotated_by_Promoter.csv`: DMP annotated promoter's information.

- `AllProbe_Related_ExpandedGenes.csv`: Expanded Genes related to all probes.

- `DMP_Related_ExpandedGenes.csv`: Expression genes related to all DMPs.

- `AllProbe_Related_Promoters.csv`: Promoters related to all probes.

- `DMP_Related_Promoters.csv`: Promoters related to all DMPs.

## Data folder:

- `AllProbe_ExpandedGene.rda`: R object for mapping result between all probes and expanded genes.

- `DMP_ExpandedGene.rda`: R object for mapping result between DMPs and expanded genes.

- `AllProbe_Promoter.rda`: R object for mapping result between all probes and promoters.

- `DMP_Promoter.rda`: R object for mapping results between DMPs and promoters.

## Figure folder:

- `AtLeast5DMPEnrichedExpandedGene.pdf`: heatmap of DNA methylation value for all genes enriched by at least 5 CpGs. Note that gene ranges are expanded to upstream 2000. 

- `CpGGeneScatterPlot_Labeled.pdf`: Scatter plot between `beta` value differential beta` value for CpGs, and log2FC from gene expression value. With all CpGs in promoter labelled.

- `AtLeast5DMPEnrichedPromoter.pdf`: Heatmap of DNA methylation value from all promoters enriched by at least 5 CpGs. Promoter range is calcuated from `promoters` function.

- `CpGGeneScatterPlot.pdf`: Scatter plot between beta value differential beta value for CpGs, and log2FC from gene expression value.
