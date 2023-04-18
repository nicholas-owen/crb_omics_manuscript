# 04 Visualisation of Differentially Methylated Positins

To visualise the DMPs we generated plots based upon the enriched data from stage 02.

- `GeneFeatureBP.R`: The function wrote to create gene `feature` plot based on a differential CpG matrix (`data.frame` object), and a `gff` (`genomation`) `GRanges` object.

Note that the promoter is calculated based on promoters function from `GenomicRanges` package, with default parameter, so the range is 2000 upstream of gene start site, and 200 downstream.

- `VolcanoPlot.R`: The script used to create `volcano` plot.

- `ExpandGeneRange.R`: A function to expand gene ranges from the GRCz11.98 annotation labelled start/end to a 2000 bp upstream promoter.

## Figure folder:

- `DMPVolcanoPlot.pdf`: Volcano plot, with top significant CpG's related mapped genes labelled.

- `GeneFeatureBarplot.pdf`: Number of DMPs mapped on corresponding gene `feature`.
