import pandas as pd
import utils as u
import matplotlib.pyplot as plt

counts = pd.read_table('../counted_reads/aggregated_unconsidering_sense.tsv')
counts_sum = counts.sum()
counts_sum = counts_sum.iloc[1:]
genes_data = counts_sum[counts_sum.index.str.startswith('Smp')]

head_data = pd.read_table('../genome_annotation/all_together_now.tsv')

plt.figure(figsize=(6, 10))
u.multibox_compare((
        head_data.gene_transcription.dropna(),
        genes_data.drop(head_data.neighbor_gene.dropna()),
        genes_data),
    ('Com Perere-3 vizinho', 'Sem Perere-3 vizinho', 'Total'), margin=3)

u.save_all_figs()
