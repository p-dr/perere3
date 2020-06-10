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
    dkey = d.stem.split('_')[0]
    print('Reading', dkey)
    counts[dkey] = pd.DataFrame()

    for c in d.glob('*.tsv'):
        counts[dkey] = counts[dkey].add(pd.read_table(c, names=['id', 'count'], index_col='id'), fill_value=0)

    counts[dkey] = counts[dkey][counts[dkey].index.str.startswith('head')]
    print(counts[dkey].shape)
    #counts[dkey] = counts[dkey][~counts[dkey].index.str.endswith('complement')]

counts = dict(sorted(counts.items()))

plt.figure(figsize=(4,10))
u.multibox_compare([c.loc[c['count'] != 0, 'count'] for c in counts.values()], labels=counts.keys())
plt.savefig('multilen_counts_compare_no_zero.png')


counts_df = pd.DataFrame()

for i_headlen in counts:
    counts_df[int(i_headlen.strip('bp'))] = counts[i_headlen].mean()
    for j_headlen in counts:
              print(i_headlen, '/', j_headlen, ':',
              counts[i_headlen].mean() / counts[j_headlen].mean())
        
print(counts_df)
counts_df.to_csv('counts_means.tsv', sep='\t')

# print('média (500bp/200bp)', (counts['500bp'] / counts['200bp']).replace([np.inf, -np.inf], np.nan).dropna().mean())
# print('(média 500bp) / (média 200bp)', counts['500bp'].mean() / counts['200bp'].mean())



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
