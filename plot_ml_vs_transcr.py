import pandas as pd
from utils import pardir, save_all_figs
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

data = pd.read_table(pardir/'genome_annotation/all_together_now.tsv')
data.dropna(inplace=True)
data = data[data.flag != 'olap']

data.plot('motherlength', 'transcription', 'scatter')
data['ml_bins'] = pd.cut(data.motherlength, [0, 750, 3150,
                                             data.motherlength.max()])
data.boxplot('transcription', 'ml_bins', showfliers=False, notch=False)

grouped = data.groupby('ml_bins')
print(grouped.corr(method="spearman"))
for group in grouped.groups:
    print(group)
exit()

save_all_figs()
plt.show()
