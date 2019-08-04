# description: Plots correlation as a function of mean and median of normalized read counts through each SRA library. Intends to show outlier effects on measured correlations as number of reads counted decreases.
# in: pardir/'genome_annotation/head_genes_correlations.tsv' pardir/'counted_reads/aggregated.tsv'
# out:
# plot:

from utils import pardir
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

LOWESS_FRAC = 1/5

correlations = pd.read_csv(pardir/('genome_annotation/head_genes_correlations' +
                           '.tsv'), sep='\t', usecols=['head_id', 'correlation'],
                           index_col='head_id')
agg_counts = pd.read_csv(pardir/'counted_reads/aggregated.tsv', sep='\t',
                         index_col='biblioteca')
agg_counts = agg_counts.loc[:, agg_counts.columns.str.startswith('head')]

# Head-gene data
hg_data = pd.concat([correlations,
                     agg_counts.mean(),
                     agg_counts.median()],
                    1, sort=False)

hg_data.dropna(inplace=True)
hg_data.columns = ['correlation', 'Mean', 'Median']

# hg_data['mean_lowess'] = lowess(hg_data.correlation, hg_data.Mean,
#                                 return_sorted=False, frac=1/3)
# hg_data['median_lowess'] = lowess(hg_data.correlation, hg_data.Median,
#                                   return_sorted=False, frac=1/3)

ax = None
for i, col in enumerate(hg_data.columns[1:]):
    ax = hg_data.plot.scatter(col, 'correlation', label=col,
                              ax=ax, color='C'+str(i), logx=True, alpha=.1)

plt.plot(*lowess(hg_data.correlation, hg_data.Mean, frac=LOWESS_FRAC).T, label='lowess mean')
plt.plot(*lowess(hg_data.correlation, hg_data.Median, frac=LOWESS_FRAC).T, label='lowess median')
plt.legend()

plt.show()

