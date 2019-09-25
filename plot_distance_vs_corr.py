# description: Plots distance from each head to its nearest or overlapped gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
from utils import pardir, save_all_figs, GFF3_COLUMNS, unfold_gff
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

#data['lowess'] = lowess(data.correlation, data.distance, .6, return_sorted=False)
head_data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv',
                          usecols=['flag', 'strand', 'distance',
                                   'correlation', 'neighbor_gene'])
head_data = head_data.dropna().reset_index()

gene_data = pd.read_table(pardir/'genome_annotation/gene_annotations.gff3',
                          names=GFF3_COLUMNS, header=None)
gene_data = unfold_gff(gene_data)

### SPLITS

## Strand
neighbor_genes = head_data.neighbor_gene.dropna()
nei_strand = gene_data.set_index('gene_id').loc[neighbor_genes].reset_index().strand
same_strand = head_data[head_data.strand == nei_strand]
diff_strand = head_data.drop(same_strand.index)

## up/down-stream
# we're restricting to same strand.
upstream = same_strand[same_strand.strand == '+'


## Overlaps or not
