# description: Plots distance from each head to its nearest or overlapped gene as a function of the transcription coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
from utils import pardir, save_all_figs, GFF3_COLUMNS, unfold_gff, show_flag
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import mannwhitneyu


head_data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')

# # Very important to drop NaN's! (we use stuff like data[data.flag != 'olap'])
# This keeps only heads with neighbor genes.
head_data = head_data.dropna().reset_index(drop=True)


# #################### SPLITS ##########################

# ##### OVERLAPS OR NOT
overlaps = head_data[(head_data.flag == 'olap') & head_data.same_strand]
noverlap = head_data.drop(overlaps.index)

# drop overlapped
head_data = head_data[head_data.flag != 'olap']

# ##### SAME/DIFF STRAND
same_strand = head_data[head_data.same_strand]
diff_strand = head_data.drop(same_strand.index)

# ##### UP/DOWN-STREAM: ----->   -->
downstream = same_strand[(same_strand.strand == '+') &
                         (same_strand.flag == 'dir')]
downstream = downstream.append(same_strand[(same_strand.strand == '-') &
                                           (same_strand.flag == 'esq')])

notdownstream = head_data.drop(downstream.index)

# ##### COMPLETE OR NOT
complete = head_data[head_data.motherlength > 3150]
notcomplete = head_data[head_data.motherlength < 750]

nolap = (complete.flag == 'olap').sum()
tot = complete.flag.count()
print(f'Quantos completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')

nolap = (notcomplete.flag == 'olap').sum()
tot = notcomplete.flag.count()
print(f'Quantos não completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')

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
         (close_promoter, not_close_promoter))

labels = ['Sobreposta ou não sobreposta ',
          'Mesma fita ou em fitas diferentes ',
          'Downstream ou não downstream ',
          'Completa ou incompleta ',
          'Promotor próximo ou promotor mais distante ']

for label, pair in zip(labels, pairs):
    a, b = [p.transcription for p in pair]
    pvalue = mannwhitneyu(a, b).pvalue
    print(label, pvalue)

    if 1:
        plt.figure()
        plt.title(label + str(pvalue))

        plt.hist(a, alpha=.5)
        plt.hist(b, alpha=.5)
        plt.legend(label[:-1].split(' ou '))
        plt.xlabel('Transcrição')
        plt.ylabel('Frequência')

if show_flag:
    plt.show()
else:
    save_all_figs()


# ################# PLOT ######################

for pair, title in zip(pairs[1:], labels[1:]):
    fig = plt.figure()
    plt.title(title)
    leg_labs = title[:-2].split(' ou ')

    for i, data in enumerate(pair):
        cor = 'C' + str(i)

        data.sort_values('distance', inplace=True)
        data.loc[:, 'lowess'] = lowess(data.transcription, data.distance, .6,
                                       return_sorted=False)

        plt.semilogx(*data[['distance', 'lowess']].values.T,
                     zorder=10, label=leg_labs[i], color=cor)

        fig = plt.semilogx(*data[['distance', 'transcription']].values.T,
                           '.', alpha=.4, color=cor)

    plt.legend()
    plt.ylabel('Transcrição')
    plt.xlabel('Distância ao gene vizinho (bp)')
    # plt.axis([1e3, 1e5, -.1, .4])

if show_flag:
    plt.show()
else:
    save_all_figs()
