import os.path

from tissue_enrichment import tscores
from tissue_enrichment import gene_alias_from_hugo


if __name__ == '__main__':
    ts = tscores.TissueScores(hugo_aliases=False)

    match, non_match = gene_alias_from_hugo.read_hugo(os.path.join('data', 'non_alt_loci_set.txt'), ts)
    idx = gene_alias_from_hugo.create_idx(match, non_match)
    gene_alias_from_hugo.save_idx_lgroup(idx, 
                                         match, 
                                         non_match, 
                                         os.path.join('data', 'idx.csv'), 
                                         os.path.join('data', 'locus_group.csv'))
