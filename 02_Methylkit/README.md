# 02 Methylkit Analysis

The second step used the Methylkit Bioconductor package (version `1.14.2`) to process bisMark mapped `BAM` files.

We applied the default function as showed in [methylkit vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) to process Bismark mapped reads into CpG information. Note that compared with default parameter, the `mincov` parameter (indicates CpG coverage) is `7`, and `minqual` (indicates quality) is `20`.

The `filterByCoverage()` is used to do filtering for CpGs with less than 7 reads, and CpGs with too much outlier number reads. The cutoff for top outlier `hi.perc` is 99.9. The beta matrix is generated without offset.

## Data Folder:

- `beta.rda`: beta matrix generated for all 6 samples, which contains __746647__ CpGs.

- `meth.rda`: meth object defined by `methylkit`, used to calcualte DMP later, which contains __746647__ CpGs.

- `myobj.rda`: raw `methylkit` object.

## Figure Folder:

- `DensityPlot.pdf`: Density of CpG methylation across all sites for all 6 samples.

- `CoveragePlot.pdf`: CpG coverage of CpG methylation across all sites for all 6 samples.

- `FilteredCoveragePlot.pdf`: CpG coverage of CpG methylation across all sites for all 6 samples after filtering and coverage adjustment.

