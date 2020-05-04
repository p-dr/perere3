# description: Plots distance from each head to its nearest or overlapping gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
import string
from utils import (pardir, save_all_figs,
                   show_flag, boxplot)
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import mannwhitneyu

plt.rcParams['font.size'] = 8

data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')
data = data[data.repetitions == 0]
# WARNING: There are missing neighbor genes. Consider these NaN values while analysing data.

COMPLETE_THRESHOLD = 3150
INCOMPLETE_THRESHOLD = 750
abc = [l + ') ' for l in string.ascii_letters]


def split(data):

    splits = dict()

    # ##### OVERLAPS OR NOT
    overlapping = data[data.relative_position == 'olap']
    # ## Heads with no neighbor gene are also not overlapping.
    # not_overlapping = data[data.relative_position.isin({'hg', 'gh'})]
    not_overlapping = data.drop(overlapping.index)
    splits[('Overlapping', 'Not overlapping')] = (overlapping, not_overlapping)

    # ##### COMPLETE OR NOT
    complete = data[data.motherlength > COMPLETE_THRESHOLD]
    incomplete = data[data.motherlength < INCOMPLETE_THRESHOLD]
    splits[('Complete', 'Incomplete')] = (complete, incomplete)

    # ## scaff.
    nolap = (complete.relative_position == 'olap').sum()
    tot = complete.relative_position.count()
    print(f'Quantos completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')

    nolap = (incomplete.relative_position == 'olap').sum()
    tot = incomplete.relative_position.count()
    print(f'Quantos não completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')
    ###

    # ##### SAME/DIFF STRAND
    same_strand = not_overlapping[not_overlapping.same_strand]
    different_strand = not_overlapping.drop(same_strand.index)
    splits[('Same strand n.o.', 'Different strand n.o.')] = (same_strand, different_strand)

    # ##### UP/DOWN-STREAM
    # (... of the neighbor gene)
    downstream = data[data.gene_stream == 'gh']
    upstream = data[data.gene_stream == 'hg']
    splits[('Downstream of NG', 'Upstream of NG')] = (downstream, upstream)

    # ##### COMPLETE OR NOT AND NOT OVERLAPPING
    complete_not_overlapping = not_overlapping[not_overlapping.motherlength > COMPLETE_THRESHOLD]
    incomplete_not_overlapping = not_overlapping[not_overlapping.motherlength < INCOMPLETE_THRESHOLD]
    splits[('Complete n.o.', 'Incomplete n.o.')] = (complete_not_overlapping, incomplete_not_overlapping)

    # ##### CLOSE PROMOTER: <--   -------->
    close_promoters = data[(data.gene_stream == 'hg') & ~data.same_strand]
    distant_promoters = data[(data.gene_stream == 'gh') & ~data.same_strand]
    splits[('Close promoters', 'Distant promoters')] = (close_promoters, distant_promoters)

    return splits


def plot_as_boxes(corr_transcr, splits, ylabel, selected=None):
    
    if selected is None:
        selected = list(range(len(splits)))

    n_splits = len(selected)

    # ## Create outer frame.
    fig, axs = plt.subplots(1, n_splits, figsize=(4 * n_splits, 4.8))
    for ax in axs:
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

    ax = fig.add_subplot(111, frameon=False)
    ax.grid(False)

    plt.tick_params(labelcolor='none',
                    top=False, bottom=False, left=False, right=False)
    print('\n***** ', ylabel, ' *****\n')
    plt.ylabel(ylabel, labelpad=20)

    # ## Plot pair comparisons.
    i = 0
    for labels, data_pair in splits.items():
        a, b = [p[corr_transcr].dropna() for p in data_pair]
        pvalue = mannwhitneyu(a, b).pvalue
        plabel = f'p-valor:\n{pvalue:.5}'
        title = ' / '.join(labels)
        median_ratio = data_pair[0][corr_transcr].median() / data_pair[1][corr_transcr].median()
        print(f'{title:<40}', '|  p-value:',
              pvalue, ['x', 'o'][pvalue < .05],
              '| Median ratio:', median_ratio, sep='\t')

        if i in selected:
            fig.add_subplot(1, len(selected), selected.index(i) + 1, frameon=False)
            plt.title(abc[selected[i]] + title)

            boxplot([a, b])

            plt.annotate(plabel, (.5, .885), xycoords='axes fraction', ha='center')
            labels = [l + '\n' + str(len(data)) for l, data in zip(labels, [a, b])]
            plt.xticks([0, 1], labels=labels)
    
        i += 1

    #plt.tight_layout()
    # plt.subplots_adjust()


def plot_as_scatter(corr_transcr, splits, ylabel):

    # ################# SCATTER PLOT ######################

    for labels, data_pair in splits.items():
        title = ' / '.join(labels)
        fig = plt.figure(figsize=(9, 4.8), dpi=200)
        plt.title(title)

        for i, data in enumerate(data_pair):
            cor = 'C' + str(i)

            data.sort_values('distance', inplace=True)
            data['lowess'] = lowess(data[corr_transcr], data.distance, .6,
                                           return_sorted=False)

            plt.semilogx(*data[['distance', 'lowess']].values.T,
                         zorder=10, label=labels[i], color=cor)

            fig = plt.semilogx(*data[['distance', corr_transcr]].values.T,
                               '.', alpha=.4, color=cor)

        plt.legend()
        plt.ylabel(ylabel)
        plt.xlabel('Distância ao gene vizinho (bp)')
        # plt.axis([1e3, 1e5, -.1, .4])


def main(corr_transcr, splits):

    if corr_transcr == 'transcription':
        ylabel = 'RPKM'
    elif corr_transcr == 'correlation':
        ylabel = 'Correlação transcricional com o gene vizinho'

    plot_as_boxes(corr_transcr, splits, ylabel)   
    #plot_as_scatter(corr_transcr, splits, ylabel)   


################# MAIN #########################
splits = split(data)

for corr_transcr in ('transcription', 'correlation'):
    main(corr_transcr, splits)

if show_flag:
    plt.show()
else:
    save_all_figs()
