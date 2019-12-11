# description: Plots distance from each head to its nearest or overlapped gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
from utils import (pardir, save_all_figs,
                   show_flag, boxplot)
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import mannwhitneyu


head_data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')

# # Very important to drop NaN's! (we use stuff like data[data.flag != 'olap'])
# This keeps only heads with neighbor genes.
head_data = head_data.dropna().reset_index(drop=True)


# ################# CORRELATION ########################
print("TABELA DE CORRELAÇÕES DE SPEARMAN")
print(head_data[['transcription',
                 'correlation',
                 'distance']].corr(method='spearman'))


# #################### SPLITS ##########################

# ##### OVERLAPS OR NOT
# overlaps = head_data[(head_data.flag == 'olap') & head_data.same_strand]
overlaps = head_data[(head_data.flag == 'olap')]
print("CORRELAÇÕES DOS OLAP")
print(overlaps[['transcription', 'correlation']].corr(method='spearman'))

# noverlap = head_data.drop(overlaps.index)
noverlap = head_data[head_data.flag != 'olap']

print(noverlap[['correlation', 'flag']])
print("CORRELAÇÕES DOS NOLAP")
print(noverlap[['transcription', 'correlation']].corr(method='spearman'))

# ##### COMPLETE OR NOT
complete = head_data[head_data.motherlength > 3150]
notcomplete = head_data[head_data.motherlength < 750]

nolap = (complete.flag == 'olap').sum()
tot = complete.flag.count()
print(f'Quantos completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')

nolap = (notcomplete.flag == 'olap').sum()
tot = notcomplete.flag.count()
print(f'Quantos não completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')


# ###################### DROP OVERLAPPED ############################ #
head_data = head_data[head_data.flag != 'olap']
# ################################################################### #


# ##### SAME/DIFF STRAND
same_strand = head_data[head_data.same_strand]
diff_strand = head_data.drop(same_strand.index)

# ##### UP/DOWN-STREAM: ----->   -->
# downstream = same_strand[(same_strand.strand == '+') &
#                          (same_strand.flag == 'dir')]
# downstream = downstream.append(same_strand[(same_strand.strand == '-') &
#                                            (same_strand.flag == 'esq')])

downstream = head_data[(head_data.strand == '+') &
                       (head_data.flag == 'dir')]
downstream = downstream.append(head_data[(head_data.strand == '-') &
                                         (head_data.flag == 'esq')])

notdownstream = head_data.drop(downstream.index)

# ##### COMPLETE OR NOT AND NOT OVERLAPPED
complete_nolap = head_data[head_data.motherlength > 3150]
notcomplete_nolap = head_data[head_data.motherlength < 750]

# ##### CLOSE PROMOTER: <--   -------->
close_promoter = diff_strand[(diff_strand.strand == '+') &
                             (diff_strand.flag == 'dir')]
close_promoter.append(diff_strand[(diff_strand.strand == '-') &
                                  (diff_strand.flag == 'esq')])

not_close_promoter = head_data.drop(close_promoter.index)


# ################# Wilcoxon #####################

pairs = ((overlaps, noverlap),
         (same_strand, diff_strand),
         (downstream, notdownstream),
         (complete, notcomplete),
         (complete_nolap, notcomplete_nolap),
         (close_promoter, not_close_promoter))

labels = ['Sobreposta ou não sobreposta ',
          'Mesma fita ou em fitas diferentes ',
          'Downstream ou upstream ',
          'Completa ou incompleta ',
          'Completa externa ou incompleta externa',
          'Promotor próximo ou promotor mais distante ']

abc = ('a) ', 'b) ', 'c) ', 'd) ')
selected = (0, 2, 3, 4)


fig, axs = plt.subplots(1, len(selected), figsize=(11, 4.8))
for ax in axs:
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

ax = fig.add_subplot(111, frameon=False)
ax.grid(False)

plt.tick_params(labelcolor='none',
                top=False, bottom=False, left=False, right=False)
plt.ylabel('Correlação transcricional com o gene vizinho', labelpad=20)

for i, label, pair in zip(range(len(labels)), labels, pairs):
    a, b = [p.correlation for p in pair]
    pvalue = mannwhitneyu(a, b).pvalue
    plabel = f'p-valor:\n{pvalue:.5}'
    print(label, pvalue)

    if i in selected:
        fig.add_subplot(1, len(selected), selected.index(i) + 1, frameon=False)
        plt.title(abc[selected.index(i)] + label)

        boxplot([a, b])

        plt.annotate(plabel, (.5, .2), xycoords='axes fraction', ha='center')
        plt.xticks([0, 1], labels=[s.capitalize()
                                   for s in label[:-1].split(' ou ')])

plt.tight_layout()
# plt.subplots_adjust()
if show_flag:
    plt.show()


# ################# PLOT ######################

for pair, title in zip(pairs[1:], labels[1:]):
    fig = plt.figure(figsize=(9, 4.8), dpi=200)
    plt.title(title)
    leg_labs = [s.capitalize() for s in title[:-1].split(' ou ')]

    for i, data in enumerate(pair):
        cor = 'C' + str(i)

        data.sort_values('distance', inplace=True)
        data.loc[:, 'lowess'] = lowess(data.correlation, data.distance, .6,
                                       return_sorted=False)

        plt.semilogx(*data[['distance', 'lowess']].values.T,
                     zorder=10, label=leg_labs[i], color=cor)

        fig = plt.semilogx(*data[['distance', 'correlation']].values.T,
                           '.', alpha=.4, color=cor)

    plt.legend()
    plt.ylabel('Correlação transcricional com o gene vizinho')
    plt.xlabel('Distância ao gene vizinho (bp)')
    # plt.axis([1e3, 1e5, -.1, .4])

if show_flag:
    plt.show()
else:
    save_all_figs()
