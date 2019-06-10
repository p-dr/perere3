from utils import read_tsv, pardir
import pandas as pd
from seaborn import heatmap
from matplotlib import pyplot as plt
import numpy as np
from plot_read_counts import plot_heatmap

ml = read_tsv(pardir/'genome_annotation/heads_motherlength.tsv', header=None,
              names=['head_id', 'ml'], index_col='head_id')

corr = read_tsv(pardir/'genome_annotation/head_genes_correlations.tsv', header=None,
                names=['head_id', 'gene_id', 'flag', 'corr'],
                index_col='head_id')

ml_corr = pd.merge(ml, corr, left_index=True, right_index=True).dropna()
#print(ml_corr)
ml_corr.plot('ml', 'corr', 'scatter', alpha=.2)
ml_corr.plot.hexbin('ml', 'corr', gridsize=15, cmap='viridis')

#======================== heatmap ========================#

# bins = 50

# ml_corr['ml_bins'] = pd.cut(ml_corr['ml'], bins)
# ml_corr['corr_bins'] = pd.cut(ml_corr['corr'], bins)
# ptab = pd.pivot_table(ml_corr, values=np.array([1]*len(ml_corr)),
#                       index='corr_bins', columns='ml_bins',
#                       aggfunc=lambda x: 1)

# ptab.fillna(0, inplace=True)
# plt.pcolor(ptab)
# plt.colorbar()
# plt.xlabel('motherlength')
# plt.ylabel('coeficiente Pearson')

ml_corr['bins'] = pd.cut(ml_corr['ml'], 3)
ml_corr.boxplot('corr', 'bins')

print(ml_corr.groupby('bins').mean())
plt.show()
