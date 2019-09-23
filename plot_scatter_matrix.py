# in:
    # pardir/'genome_annotation/head_genes_correlations_unconsidering_sense.tsv'
    # pardir/'genome_annotation/head_annotations.gff3'
    # pardir/'counted_reads/aggregated_unconsidering_sense.tsv

import pandas as pd
import matplotlib.pyplot as plt
from utils import pardir, save_all_figs, GFF3_COLUMNS, unfold_gff, plot_box, show_flag

corr_data = pd.read_table(pardir/'genome_annotation/head_genes_correlations_unconsidering_sense.tsv')
corr_data.set_index('head_id', inplace=True)

gff_data = pd.read_table(pardir/'genome_annotation/head_annotations.gff3', header=None,
                         names=GFF3_COLUMNS)
gff_data = unfold_gff(gff_data)
gff_data.set_index('gene_id', inplace=True)
gff_data.motherlength = gff_data.motherlength.astype(dtype='int')

count_data = pd.read_table(pardir/'counted_reads/aggregated_unconsidering_sense.tsv')
count_data = count_data.loc[:, count_data.columns.str.startswith('head')]
count_data = count_data.sum()

data = pd.concat([gff_data, corr_data, count_data], axis=1)

# PLOT ALL!
pd.plotting.scatter_matrix(data.drop(['end', 'start'], 1))
data.plot.scatter('start', 0, alpha=.2)
data.plot.scatter('distance', 0, alpha=.2)
plot_box(data[data.distance > 0], 'distance', 0, bins=50)

genome_map = data.sort_values('start')[['start', 0]]
genome_map['start'] = pd.cut(genome_map.start, 1e4)
genome_map = genome_map.groupby('start').sum()
l = len(genome_map[0])
plt.figure()
plt.pcolor(genome_map[0].values.reshape(int(l**.5), int(l**.5)))
plt.colorbar()
plt.figure()

plt.pcolor(data.corr(method='spearman'))
plt.figure()
data[0].hist()

if not show_flag:
    save_all_figs()
plt.show()
