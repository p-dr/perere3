# description: Plots heads' to nearest gene(s) expression correlation as a function of heads' motherlengths.
# in: pardir/'genome_annotation/heads_motherlength.tsv' pardir/'genome_annotation/head_genes_correlations.tsv'
# out: 
# plot: 

from utils import read_tsv, pardir, save_all_figs
import pandas as pd
from seaborn import heatmap
from matplotlib import pyplot as plt
import numpy as np
from plot_read_counts import plot_heatmap

ml = read_tsv(pardir/'genome_annotation/heads_motherlength.tsv', header=None,
              names=['head_id', 'ml'], index_col='head_id')

corr = read_tsv(pardir/'genome_annotation/head_genes_correlations.tsv',
                usecols=['head_id', 'correlation'],
                index_col='head_id')

ml_corr = pd.merge(ml, corr, left_index=True, right_index=True).dropna()
#print(ml_corr)
ml_corr.plot('ml', 'correlation', 'scatter', alpha=.2)
ml_corr.plot.hexbin('ml', 'correlation', gridsize=15, cmap='viridis')

#======================== heatmap ========================#

# bins = 50

# ml_corr['ml_bins'] = pd.cut(ml_corr['ml'], bins)
# ml_corr['corr_bins'] = pd.cut(ml_corr['correlation'], bins)
# ptab = pd.pivot_table(ml_corr, values=np.array([1]*len(ml_corr)),
#                       index='corr_bins', columns='ml_bins',
#                       aggfunc=lambda x: 1)

# ptab.fillna(0, inplace=True)
# plt.pcolor(ptab)
# plt.colorbar()
# plt.xlabel('motherlength')
# plt.ylabel('coeficiente Pearson')

ml_corr['bins'] = pd.cut(ml_corr['ml'], 15)
ml_corr.boxplot('correlation', 'bins')

print(ml_corr.groupby('bins').mean())
save_all_figs()
plt.show()
