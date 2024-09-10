# SOF_mutants_IkBa_intestinal_stemness

This repository includes scripts required for the bulk RNAseq and ChIPseq IkBa data analysis included in Alvarez-Villanueva et al. (Under peer-review). All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Scripts from 1 to 4.
- Downstream analysis: to conduct differential expression analysis and functional analysis (GSEA). Scripts 5 and 6.

To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE206515 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206515]. All required scripts were executed using Singularity images (v3.8.3) per required tool.

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded. 
NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 106, hg38).

## ChIPseq data analysis folder

XXXXX
