import os.path

import tissue_enrichment.tissue_scores_from_gct


if __name__ == '__main__':
    import os.path
    raw_path = os.path.join('data', 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
    save_path = os.path.join('data', 'gene_median_tpm.csv')
    tissue_enrichment.tissue_scores_from_gct.process_raw(raw_path, save_path)
