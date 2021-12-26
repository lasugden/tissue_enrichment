# tissue_enrichment

Given an  input set of gene names and their associated weights, such as from the output of SWIFR or another mechanism for identifying selective sweeps, identify the tissue with the greatest enhancement, along with an associated p-value via bootstrapping.

Sugden Lab, Duquesne University, 2022

### Steps for use:
1. Download ```GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct``` from https://gtexportal.org/home/datasets
2. Convert raw data into an appropriate CSV format for reading
3. Run ```tissue_weighter.py``` to identify tissue weights.