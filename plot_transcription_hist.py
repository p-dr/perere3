import pandas as pd
from utils import pardir, save_all_figs
import matplotlib.pyplot as plt

inpath = pardir/'counted_reads/aggregated_unconsidering_sense.tsv'
data = pd.read_table(inpath).sum()[1:]
heads = data.loc[data.index.str.startswith('head')]
genes = data.drop(heads.index)
data = pd.concat([heads.astype(float), genes.astype(float)], 1)
data.columns = ['heads', 'genes']
# zeros = int((data < 5).sum())
# print(data)
# l = int(len(data))

# print(f'\n\nporcentagem de zeros: {zeros} / {l} = {100 * zeros / l}%')
# plt.hist(data.genes, bins=50, alpha=.2)
# plt.hist(data.heads, bins=50)
heads.hist(bins=30, log=True)
plt.title('Distribuição de coeficientes de transcrição entre as cópias')
plt.xlabel('Coeficiente de transcrição')
plt.ylabel('Frequência')

save_all_figs()
plt.show()
