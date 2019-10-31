import pandas as pd
from utils import pardir, save_all_figs, show_flag
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')
data.dropna(inplace=True)
data = data[data.flag != 'olap']

data.plot('motherlength', 'transcription', 'scatter', logy=True, alpha=.2)
data['ml_bins'] = pd.cut(data.motherlength, [0, 750, 3150,
                                             data.motherlength.max()])
boxprops = dict(color='black')
medianprops = dict(color='C1')
data.boxplot('transcription', 'ml_bins', showfliers=False,
             boxprops=boxprops, medianprops=medianprops,
             widths=.75, figsize=(3.5, 4.8))
plt.xlabel('Comprimento da região invariante (pb)')
plt.ylabel('RPKM')

grouped = data.groupby('ml_bins')
print('p-value:', mannwhitneyu(grouped.groups[pd.Interval(0, 750, closed='right')],
                               grouped.groups[pd.Interval(3150, 3193, closed='right')]))

plt.title('Transcrição em função do com-\nprimento da região invariante')
plt.suptitle('')

if show_flag:
    plt.show()
else:
    save_all_figs()
