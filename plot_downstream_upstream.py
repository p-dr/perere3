# description: Plots distance from each head to its nearest or overlapped gene as a function of the correlation coeficient between their normalized read counts.
# in: pardir/'genome_annotation/all_together_now.tsv'
# out: 

import pandas as pd
from utils import pardir, save_all_figs, show_flag
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


head_data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')

# # Very important to drop NaN's! (we use stuff like data[data.flag != 'olap'])
# This keeps only heads with neighbor genes.
head_data = head_data.dropna().reset_index(drop=True)


# ################# CORRELATION ########################
print("TABELA DE CORRELAÇÕES DE SPEARMAN")
print(head_data[['transcription', 'correlation', 'distance']].corr(method='spearman'))


# #################### SPLITS ##########################

# drop overlapped
head_data = head_data[head_data.flag != 'olap']

# ##### SAME/DIFF STRAND
same_strand = head_data[head_data.same_strand]
diff_strand = head_data.drop(same_strand.index)

# ##### UP/DOWN-STREAM: ----->   -->
# downstream = same_strand[(same_strand.strand == '+') &
#                          (same_strand.flag == 'dir')]
# downstream = downstream.append(same_strand[(same_strand.strand == '-') &
#                                            (same_strand.flag == 'esq')])

# upstream = same_strand[(same_strand.strand == '+') &
#                          (same_strand.flag == 'esq')]
# upstream = upstream.append(same_strand[(same_strand.strand == '-') &
#                                            (same_strand.flag == 'dir')])

# ##==============================================================
downstream = head_data[(head_data.strand == '+') &
                       (head_data.flag == 'dir')]
downstream = downstream.append(head_data[(head_data.strand == '-') &
                                         (head_data.flag == 'esq')])

upstream = head_data.drop(downstream.index)



# ################# Wilcoxon #####################

thresholds = (1e3, 1e4, 2e4, upstream.append(downstream).distance.max())

abc = ('a) ', 'b) ', 'c) ')

fig, axs = plt.subplots(1, len(thresholds), figsize=(11, 4.8))
for ax in axs:
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none',
                top=False, bottom=False, left=False, right=False)
plt.ylabel('Correlação transcricional com o gene vizinho')

for i, thresh in enumerate(thresholds):
    # Data selection
    a, b = [p.loc[p.distance <= thresh].correlation
            for p in [downstream, upstream]]

    pvalue = mannwhitneyu(a, b).pvalue
    plabel = f'p-valor:\n{pvalue:.5f}'
    label = f"Distância < {thresh:.0f}"

    print(label, pvalue, 'Medians:', a.median(), b.median())

    fig.add_subplot(1, len(thresholds), i + 1, frameon=False)
    plt.title(label)

    plt.boxplot([a, b], widths=.75)
    plt.annotate(plabel, (.5, .2), xycoords='axes fraction', ha='center')
    plt.xticks([1, 2], labels=['Downstream', 'Upstream'])

plt.tight_layout()


# ======================= P VS. DIS ==========================

plt.figure(figsize=(11, 4.8))
plt.subplot(211)

xdistances = range(100, int(upstream.distance.max()), 100)
ypvalues = [mannwhitneyu(downstream.loc[downstream.distance <= thresh].correlation,
                         upstream.loc[upstream.distance <= thresh].correlation).pvalue
            for thresh in xdistances]

updataamounts = [len(upstream.loc[upstream.distance <= thresh]) for thresh in xdistances]
downdataamounts = [len(downstream.loc[downstream.distance <= thresh]) for thresh in xdistances]

plt.plot(xdistances, ypvalues)
plt.semilogx()
plt.ylabel('p-valor')

# ----------------------------------
plt.subplot(212)
plt.plot(xdistances, downdataamounts, label="Downstream")
plt.plot(xdistances, updataamounts, label="Upstream")
plt.legend()

plt.semilogx()
plt.xlabel('Distância ao vizinho (pb)')
plt.ylabel('Quantidade de pontos')

# ===========================================
# Só na faixa selecionada
plt.figure(dpi=200)
limits = 1e3, 1e4
print(i.distance.between(*limits))
plt.boxplot([[i.correlation.loc[i.distance.between(*limits)]]
             for i in (upstream, downstream)])

# ===========================================
if show_flag:
    plt.show()
else:
    save_all_figs()
