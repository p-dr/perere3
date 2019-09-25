# description: Plots distance from each head to its nearest or overlapped gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
from utils import pardir, save_all_figs, GFF3_COLUMNS, unfold_gff
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

#data['lowess'] = lowess(data.correlation, data.distance, .6, return_sorted=False)
head_data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv',
                          usecols=['flag', 'strand', 'distance', 'correlation'])
gene_data = pd.read_table(pardir/'genome_annotation/gene_annotations.gff3',
                          names=GFF3_COLUMNS, header=None)
gene_data = unfold_gff(gene_data)
print(gene_data)


### SPLITS

## Strand
ax = data.plot('distance', 'correlation', 'scatter', alpha=.2)
data.sort_values('distance').plot.line('distance', 'lowess', ax=ax, color='C1')
plt.show()
