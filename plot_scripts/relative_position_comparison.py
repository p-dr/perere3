# description: Compare p-value between different head-gene configurations, considering strands.
# in: pardir/'genome_annotation/all_together_now.tsv'

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils import pardir, show_flag, box_compare, multibox_compare, save_all_figs, get_subsets
import utils as u

d = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')
# d = d[d.repetitions == 0]


def compare(subsets, corr_transcr):
    plt.figure(figsize=(4, 8))
    plt.title(corr_transcr)
    plt.ylabel(corr_transcr)
    pp = multibox_compare(*zip(*[(df.dropna(), key) for key, df in subsets.items()]))
    for labels, p in pp.items():
        seta, setb = labels
        print(labels, subsets[seta].describe(), subsets[setb].describe(), p, sep='\n\t')


def main():
    #print('Quantidades das populações')
    #[print(name, s) for name, s in u.get_subsets(d).items()]

    for corr_transcr in ('transcription', 'complement_transcription',
                         'gene_transcription', 'gene_complement_transcription',
                         'correlation', 'complement_correlation',
                         'motherlength', 'repetitions'):
        u.print_header(corr_transcr)
        # D = d[d[corr_transcr] != 0].dropna(subset=[corr_transcr])  # drop zeros
        D = d.dropna(subset=[corr_transcr])
        # print(f'Comprimento total: {d.shape[0]} | Comprimento coluna: {D.shape[0]}')

        if 'gene' in corr_transcr:
            D = D.sort_values('distance')  # Keep only closest copy to each G
            D = D.drop_duplicates(subset='neighbor_gene')

        subsets = get_subsets(D, corr_transcr)
        compare(subsets, corr_transcr)

    if show_flag:
        plt.show()
    else:
        save_all_figs()


if __name__ == '__main__':
    main()
