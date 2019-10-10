import pandas as pd
import matplotlib.pyplot as plt
from utils import pardir
import numpy as np


# def transc_map(data):
#     seg_data = pd.cut(data

inpath = pardir/'genome_annotation/all_together_now.tsv'

data = pd.read_table(inpath)
groups = data.groupby('seqid')

for name, group in groups:
    nonzero = group.transcrição.astype('bool').sum()
    if nonzero < 10:
        continue
    group = group.sort_values('start')
    print('nonzero:', nonzero)
    print(group[group.transcrição > 0])
    group.transcrição[group.transcrição > 0].hist(log=True)
    plt.figure()
    plt.bar(group.start.values, group.transcrição.values, width=1000,
            log=True)
    plt.title(name)
    plt.show()
