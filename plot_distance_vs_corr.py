# description: Plots distance from each head to its nearest or overlapped gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/head_genes_correlations.tsv'
# out: 

from pandas import read_csv
from utils import pardir
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

data = read_csv(pardir/'genome_annotation/head_genes_correlations.tsv', usecols=['distance', 'correlation'], sep='\t')
data['lowess'] = lowess(data.correlation, data.distance, .6, return_sorted=False)

ax = data.plot('distance', 'correlation', 'scatter', alpha=.2)
data.sort_values('distance').plot.line('distance', 'lowess', ax=ax, color='C1')
plt.show()
