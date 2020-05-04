# description: Compares genes' to gene complements' transcription in each reads library.
# in: pardir/'counted_reads/aggregated.tsv'
# plot: 
import sys
import utils as u
sys.path.append(str(u.scripts_dir))
from correlate_heads_to_near_genes import out_rpkm_by_lib as rpkm_by_lib

import pandas as pd
import matplotlib.pyplot as plt

d = pd.read_table(rpkm_by_lib, index_col='library')
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
