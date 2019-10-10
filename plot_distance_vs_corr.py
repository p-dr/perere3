# description: Plots distance from each head to its nearest or overlapped gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
from utils import pardir, save_all_figs, GFF3_COLUMNS, unfold_gff, show_flag
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import mannwhitneyu


head_data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv',
                          usecols=['flag', 'strand', 'distance',
                                   'correlation', 'neighbor_gene',
                                   'motherlength'])
# # Very important to drop NaN's! (we use stuff like data[data.flag != 'olap'])
head_data = head_data.dropna().reset_index()

gene_data = pd.read_table(pardir/'genome_annotation/gene_annotations.gff3',
                          names=GFF3_COLUMNS, header=None)
gene_data = unfold_gff(gene_data)

# ## SPLITS

# # Strand
neighbor_genes = head_data.neighbor_gene.dropna()
nei_strand = gene_data.set_index('gene_id').loc[neighbor_genes].reset_index().strand
same_strand = head_data[head_data.strand == nei_strand]
diff_strand = head_data.drop(same_strand.index)

# # up/down-stream
# we're restricting to same strand.
downstream = same_strand[(same_strand.strand == '+') &
                       (same_strand.flag == 'dir')]
downstream = downstream.append(same_strand[(same_strand.strand == '-') &
                                       (same_strand.flag == 'esq')])
# downstream = same_strand[(same_strand.strand == '+') &
#                        (same_strand.flag == 'dir')]
# downstream = downstream.append(same_strand[(same_strand.strand == '-') &
#                                        (same_strand.flag == 'esq')])
# 
upstream = same_strand.drop(downstream.index)[same_strand.flag != 'olap']
upstream.dropna(inplace=True)

# # Overlaps or not
overlaps = same_strand[same_strand.flag == 'olap']
noverlap = same_strand.drop(overlaps.index)
# ((overlaps, noverlap), 'Ovelaps/Does not overlap')

# # Complete or not
complete = same_strand[same_strand.motherlength > 3150]
notcomplete = same_strand[same_strand.motherlength < 750]

nolap = (complete.flag == 'olap').sum()
tot = complete.flag.count()
print(f'Quantos completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')

nolap = (notcomplete.flag == 'olap').sum()
tot = notcomplete.flag.count()
print(f'Quantos não completos sobrepõem: {nolap} / {tot} = {100 * nolap / tot}%')

    # Drop overlapped
complete = complete.loc[complete.flag != 'olap']
notcomplete = notcomplete.loc[notcomplete.flag != 'olap']

# ### Wilcoxon
pairs = ((overlaps, noverlap),
         (same_strand, diff_strand),
         (upstream, downstream),
         (complete, notcomplete))

labels = ['Sobreposta ou não: ',
          'Mesma fita ou não: ',
          'Jusante ou motante: ',
          'Completa ou não: ']

for label, pair in zip(labels, pairs):
    a, b = [p.correlation for p in pair]
    pvalue = mannwhitneyu(a, b).pvalue
    print(label, pvalue)

    plt.figure()
    plt.title(label + str(pvalue))
    plt.hist(a, alpha=.5)
    plt.hist(b, alpha=.5)

save_all_figs()
plt.show()


# ## PLOT
for pair, title in [(((same_strand, 'Same'), (diff_strand, 'Different')),
                     'Same/different strand'),
                    (((upstream, 'Upstream'), (downstream, 'Downstream')),
                     'Upstream vs downstream'),
                    (((complete, 'Completas'), (notcomplete, 'Não completas')),
                     'Completas vs incompletas')]:

    fig = plt.figure()
    plt.title(title)
    first = True

    for data, label in pair:
        data = data[data.flag != 'olap']
        data.sort_values('distance', inplace=True)
        data.loc[:, 'lowess'] = lowess(data.correlation, data.distance, .6,
                                return_sorted=False)
        cor = plt.semilogx(*data[['distance', 'lowess']].values.T,
                           label=label, zorder=10)[0]._color
        if first:
            ax = plt.axis()
            first = False
        fig = plt.semilogx(*data[['distance', 'correlation']].values.T,
                           '.', alpha=.4, color=cor)

    plt.axis([1e3, 1e5, -.1, .4])
    plt.legend()

    # pair_df = pd.concat([data.reset_index().correlation for data, label in pair], 1)
    # print(pair_df)
    # pair_df.hist(label=title)

if not show_flag:
    save_all_figs()
plt.show()
