import pandas as pd
import utils as u
import matplotlib.pyplot as plt

counts = pd.read_table('../counted_reads_bak/aggregated_unconsidering_sense.tsv')
head_data = pd.read_table('../genome_annotation/all_together_now.tsv')


def main(func):
    counts_sum = func(counts)
    counts_sum = counts_sum.iloc[1:]
    genes_data = counts_sum[counts_sum.index.str.startswith('Smp')]
    genes_data = genes_data.iloc[~genes_data.index.str.endswith('complement')]

    with_perere = head_data.gene_transcription.dropna().drop_duplicates()
    lone_genes = genes_data.drop(head_data.neighbor_gene.dropna().unique()).dropna()

    # print(
    #    with_perere,
    #    lone_genes,
    #    genes_data
    # )
       
    print(head_data.neighbor_gene.duplicated().sum())
    plt.figure(figsize=(6, 10))
    u.multibox_compare((
            with_perere,
            lone_genes,
            genes_data),
        ('Com Perere-3 vizinho', 'Sem Perere-3 vizinho', 'Total'), margin=3)


if __name__ ==  '__main__':
    for func in (pd.DataFrame.sum, pd.DataFrame.max):
        main(counts)
    u.save_all_figs()
