import pandas as pd
import matplotlib.pyplot as plt
from utils import show_flag, save_all_figs

a = pd.read_table('../alinhamentos/heads_vs_heads.bl', header=None)
lenghts = a[9][a[9] < 1000]
lenghts.hist(bins=500)
plt.title('Distribuição dos comprimentos entre alinhamentos sonda-sonda')
plt.xlabel('Comprimento do alinhamento sonda-sonda (pb)')
plt.ylabel('Quantidade de alinhamentos')

if show_flag:
    plt.show()

save_all_figs()
