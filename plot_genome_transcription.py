import pandas as pd
import matplotlib.pyplot as plt
from utils import pardir


def transc_map(data):
    seg_data = pd.cut(data

inpath = pardir/'genome_annotation/all_together_now.tsv'

data = pd.read_table(inpath)
groups = data.groupby('seqid')

for i, group in groups:
    print(group)
