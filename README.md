# tissue_enrichment

Given an  input set of gene names and their associated weights, such as from the output of SWIFR or another mechanism for identifying selective sweeps, identify the tissue with the greatest enhancement, along with an associated p-value via bootstrapping.

Sugden Lab, Duquesne University, 2022


### Steps for use:
1. Download ```GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct``` from https://gtexportal.org/home/datasets

1. Download ```Total Approved Symbols``` under ```Statistics``` in TXT format (the blue button in the right column at the bottom of the first block) from https://www.genenames.org/download/statistics-and-files/

1. Convert raw data into an appropriate CSV format for reading. Run script ```1_process_raw_gct_to_graph.py```.

1. Check that different tissues are relatively equally distributed. See conclusion 1 below. Run script ```2_check_tissue_distances.py```

1. Create aliases for gene names from HUGO. Run script ```3_create_hugo_aliases.py```


### Library files
```gene_alias_from_hugo.py```: Read in Hugo data, combine aliases, and create an index that can be searched to identify the best matching ensembl_gene_id.

```stats.py```: UNFINISHED-- primary output file for website

```tissue_scores_from_gct.py```: Code to read GTEX file to get gene transcript levels, convert to a graph, and simplify titles for future processing. Also standardizes edges.

```tscores.py```: Class to read in GTEX edge weights as well as hugo files for searching. Also includes code for randomization of statistics.


### Pre-Analysis conclusions:
1. Tissues show few clusters that might cause problems with ranking, with the exception perhaps of skin exposed and not exposed to the sun. Therefore, those two tissues are combined in ```tissue_scores_from_gct.py``` for all future analyses. ![image info](./graphs/tissue_umap.png)

1. Testis shows overexpression and a massive number of unique genes (found from compare_tissue_gene_expression.py)