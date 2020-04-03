# in: pardir/'counted_reads/aggregated.tsv'
# plot: 

import utils as u
import pandas as pd
import matplotlib.pyplot as plt
from  correlate_heads_to_near_genes import out_aggregated_counts

d = pd.read_table(out_aggregated_counts, index_col='biblioteca')
d = d.loc[:, d.columns.str.contains('Smp')]  # Restrict to genes.

comp = d.loc[:, d.columns.str.endswith('_complement')]
direct = d.drop(comp.columns, axis=1)
comp = comp.T
direct = direct.T

plt.figure(figsize=(comp.shape[1]*3, 5))

for i, lib in enumerate(comp.columns):
    plt.subplot(1, comp.shape[1], i+1)
    plt.title(lib)
    plt.ylabel('RPKM')
    u.multibox_compare([direct[lib], comp[lib]],
                       labels=('Genes', 'Complementares'))
plt.tight_layout()
u.save_all_figs()
