import utils as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

count_flags = [
    '__no_feature',
    '__ambiguous',
    '__too_low_aQual',
    '__not_aligned',
    '__alignment_not_unique',
]

dirs = u.pardir.glob('[0-9]*_counted_reads')

counts = {}

for d in dirs:
    dkey = int(d.stem.split('bp_')[0])
    print('Reading', dkey)
    counts[dkey] = pd.DataFrame()

    for c in d.glob('*.tsv'):
        counts[dkey] = counts[dkey].add(pd.read_table(c, names=['id', 'count'], index_col='id'), fill_value=0)

    counts[dkey] = counts[dkey][counts[dkey].index.str.startswith('head')]
    #counts[dkey] = counts[dkey][~counts[dkey].index.str.endswith('complement')]

counts = dict(sorted(counts.items()))

if u.args.plot:
    plt.figure(figsize=(5,10))
    u.multibox_compare([c.loc[c['count'] != 0, 'count'] for c in counts.values()], labels=counts.keys())
    plt.savefig('multiheadlen/counts.png')
    plt.show()

counts_df = pd.DataFrame()

for headlen, count in counts.items():
    count = count[count['count'] != 0]
    counts_df[headlen] = count['count'].describe()
    # for j_headlen in counts:
    #           print(i_headlen, '/', j_headlen, ':',
    #           counts[i_headlen].mean() / counts[j_headlen].mean())

counts_df = counts_df.T
counts_df.index.name = 'headlen'
print(counts_df)
counts_df.to_csv('counts_stats.tsv', sep='\t')

# # counts = {}
# # 
# # for dir in dirs:
# #     u.print_header(dir)
# #     counts[dir.stem] = {}
# # 
# #     for c in dir.glob('*.tsv'):
# #         print(c)
# #         counts[dir.stem].update({c.stem: pd.read_table(c)})
# 
# counts = {
#     dir.stem: {
#         c.stem: pd.read_table(c, names=['id', 'count'], index_col='id')
#         for c in dir.glob('*.tsv')
#     }
#     for dir in dirs
# }
# 
# for table_name in counts[dirs.pop().stem]:  # pop is too agressive :c
#     print(table_name)
#     print((
#         counts['500bp_counted_reads'][table_name] /
#         counts['200bp_counted_reads'][table_name]).replace((np.inf, -np.inf), np.nan).dropna().mean())
