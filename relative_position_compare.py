# description: compare p-value between different head-gene configurations, considering strands.
# in: pardir/'genome_annotation/all_together_now.tsv'

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils import pardir, show_flag, box_compare, save_all_figs, get_subsets

d = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')
d = d[d.repetitions == 0]


def main(subsets):
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
                pvalue = box_compare(subset_i.dropna(),
                                     subset_j.dropna(),
                                     labels=(label_i, label_j))
            else:
                pvalue = 0

            line.append(pvalue)
            n += 1
        pmatrix.append(line)

    plt.tight_layout()
    
    ### PLOT MEDIANS COMPARISON MATRIX
    plt.figure()
    plt.title('log da razão das medianas (positivo => coluna maior)')
        
    medians = np.array([subset.median() for subset in subsets.values()])
    print(medians)
    median_matrix = np.log(medians/ medians[:, np.newaxis])
    median_matrix = pd.DataFrame(median_matrix,
                                 columns=subsets.keys(),
                                 index=subsets.keys())
    sns.heatmap(median_matrix, annot=True)
    
    ### PLOT P-VALUES MATRIX
    plt.figure()
    plt.title('p-valores')
    pmatrix = pd.DataFrame(pmatrix,
                           columns=subsets.keys(),
                           index=subsets.keys())
    sns.heatmap(pmatrix, annot=True)


if __name__ == '__main__':
    for corr_transcr in ('transcription', 'correlation'):
        print('\n', corr_transcr.upper(), '\n')
        subsets = get_subsets(d, corr_transcr)
        main(subsets)

    if show_flag:
        plt.show()
    else:
        save_all_figs()
