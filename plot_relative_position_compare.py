# description: Compare p-value between different head-gene configurations, considering strands.
# in: pardir/'genome_annotation/all_together_now.tsv'

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils import pardir, show_flag, box_compare, multibox_compare, save_all_figs, get_subsets
import utils as u

d = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')
d = d[d.repetitions == 0]


def main(subsets, corr_transcr):
    plt.figure(figsize=(12, 12))
    pmatrix = []
    tri_ind = np.triu(np.arange(16).reshape(4, 4), k=1)
    tri_ind = tri_ind[tri_ind != 0]
    print(tri_ind)

    n = 0
    for label_i, subset_i in subsets.items():
        line = []

        for label_j, subset_j in subsets.items():

            if n in tri_ind:
                plt.subplot(4, 4, n+1)
                print('lengths:', len(subset_i), len(subset_j))
                pvalue = box_compare(subset_i.dropna(),
                                     subset_j.dropna(),
                                     labels=(label_i, label_j))
            else:
                pvalue = 0

            line.append(pvalue)
            n += 1
        pmatrix.append(line)

    plt.tight_layout()

    ### PLOT BOXES TOGETHER
    plt.figure(figsize=(4, 8))
    plt.title(corr_transcr)
    plt.ylabel(corr_transcr)
    #multibox_compare(*list(zip(*[(df.replace(0, pd.np.nan).dropna(), key) for key, df in subsets.items()])))
    multibox_compare(*list(zip(*[(df.dropna(), key) for key, df in subsets.items()])))

    ### PLOT MEDIANS COMPARISON MATRIX
    # plt.figure()
    # plt.title('log da razão das medianas (positivo => coluna maior)')
    #     
    # medians = np.array([subset.median() for subset in subsets.values()])
    # print(medians)
    # median_matrix = np.log(medians/ medians[:, np.newaxis])
    # median_matrix = pd.DataFrame(median_matrix,
    #                              columns=subsets.keys(),
    #                              index=subsets.keys())
    # sns.heatmap(median_matrix, annot=True)

    ### PLOT MEANS COMPARISON MATRIX
    plt.figure()
    plt.title('log da razão das meanas (positivo => coluna maior)')

    means = np.array([subset.mean() for subset in subsets.values()])
    print(means)
    mean_matrix = np.log(means/ means[:, np.newaxis])
    mean_matrix = pd.DataFrame(mean_matrix,
                                 columns=subsets.keys(),
                                 index=subsets.keys())
    sns.heatmap(mean_matrix, cmap='RdYlBu', annot=True)

    ### PLOT P-VALUES MATRIX
    plt.figure()
    plt.title('p-valores')
    pmatrix = pd.DataFrame(pmatrix,
                           columns=subsets.keys(),
                           index=subsets.keys())
    sns.heatmap(pmatrix, annot=True)


if __name__ == '__main__':
    print('Quantidades das populações')
    [print(name, s) for name, s in u.get_subsets(d).items()]

    for corr_transcr in ('transcription', 'complement_transcription',
                         'gene_transcription', 'gene_complement_transcription',
                         'correlation', 'complement_correlation'):
        u.print_header(corr_transcr)
        D = d[d[corr_transcr] != 0].dropna(subset=[corr_transcr])
        print(f'Retirada de zeros\nComprimento antes: {d.shape[0]} | Depois: {D.shape[0]}')

        if 'gene' in corr_transcr:
            D = D.sort_values('distance')  # Keep only closest copy to each G
            D = D.drop_duplicates(subset='neighbor_gene')

        subsets = get_subsets(D, corr_transcr)
        main(subsets, corr_transcr)

    if show_flag:
        plt.show()
    else:
        save_all_figs()
