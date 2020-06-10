import utils as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dirs = u.pardir.glob('[0-9]*_counted_reads')
counts = {}

for d in dirs:
    dkey = int(d.stem.split('bp_')[0])
    print('Reading', dkey)
    counts[dkey] = pd.read_table(d/'parsed_data/heads_rpkm.tsv')

counts = dict(sorted(counts.items()))

for transcr in ('transcription', 'complement_transcription'):
    df = pd.DataFrame()
    u.print_header(transcr)

    for headlen, rpkm in counts.items():
        rpkm = rpkm[transcr][rpkm[transcr] != 0]
        df[headlen] = rpkm.describe()
    df.columns.name = 'headlen'
    # df.index.name = 'stat'
    df = df.T
    print(df)
    df.to_csv(f'multiheadlen_{transcr}.tsv', sep='\t')

    if u.args.plot:
        plt.figure(figsize=(4,10))
        # u.multibox_compare([c[transcr][c[transcr] != 0] for c in counts.values()], labels=counts.keys())
        u.multibox_compare([c[transcr] for c in counts.values()], labels=counts.keys())
        plt.savefig(f'multiheadlen/{transcr}.png')
        plt.show()

