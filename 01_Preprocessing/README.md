# 01 Preprocessing RRBS data

This folder contains codes related to data preprocessing, everything started with the folder indicated through the variable `FileFolder` in `main.R`.

Preprocessing of the data will utilise: `bcl2fastq`, `trim_galore`, `parallel`, `bismark`, and  `samtools`.

Apart from above software, two NuGEN company specific scripts are used. The whole analysis process is followed by [NuGEN company provided pipeline](https://github.com/nugentechnologies/NuMetRRBS).

`./strip_bismark_sam.sh` was obtained from https://github.com/nugentechnologies/NuMetRRBS/blob/master/strip_bismark_sam.sh

`./nudup.py` was obtained from https://github.com/tecangenomics/nudup/blob/master/nudup.py

`./trimRRBSdiversityAdaptCustomers.py` was obtained from https://github.com/nugentechnologies/NuMetRRBS/blob/master/trimRRBSdiversityAdaptCustomers.py

The zebrafish genome was downloaded from Ensembl, [version 98](ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz)
