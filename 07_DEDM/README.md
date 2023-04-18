# 07 Differential Expression and Differential Methylation analysis

This is the scripts to find both differential methylation and differential expression genes.

The cutoff for Differential Expression Gene is: `padj: 0.05, |log2FC| >= 2, baseMean >=10`

The cutoff for methylation is: `q value 0.05, |meth.diff| >= 15`

- `PlotDEDMHeatmap.R`: A script to plot heatmap, it requires joined DEDM object, and TPM matrix users could specify labelName parameter to just show some gene's name instead of all.
Note: In lines `63` and `55` in main.R script, please reset the genes any _gene of interest_ to plot. Please ensure the gene you want to labeled can be find in `./csv/DEDM.csv`.

## CSV folder:

- `DEDM.csv`: Information for differential methylation and differential expression genes.

- `StrictDEDM.csv`: Save as above, but only genes with at least 5 CpGs enriched.


## Figure folder:

- `DEDMHeatmap_custom_label.pdf`: Heatmap for DEDM csv file

- `StrictDEDMHeatmap_custom_label.pdf`: Heatmap for StrictDEDM csv file

- `HyperVennPlot.pdf`: Venn plot between Hyper genes and Up- and Down- regulated.

- `StrictDEDMHeatmap_full_label.pdf`: Heatmap with all genes labelled

- `HypoVennPlot.pdf`: Venn plot between Hypo genes and Up- and Down- regulated. 

