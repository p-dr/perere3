import pandas as pd
from utils import pardir, save_all_figs
import matplotlib.pyplot as plt

inpath = pardir/'genome_annotation/all_together_now.tsv'
data = pd.read_table(inpath, usecols=['transcrição'])

zeros = int((data < 5).sum())
print(data)
l = int(len(data))

print(f'\n\nporcentagem de zeros: {zeros} / {l} = {100 * zeros / l}%')
data.hist(bins=20, log=True)
save_all_figs()
plt.show()
