import utils as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dirs = u.pardir.glob('[0-9]*_counted_reads')
counts = {}

for d in dirs:
    dkey = d.stem.split('_')[0]
    print('Reading', dkey)
    counts[dkey] = pd.read_table(d/'parsed_data/heads_rpkm.tsv')

counts = dict(sorted(counts.items()))
plt.figure(figsize=(4,10))

u.multibox_compare([c.transcription[c.transcription != 0] for c in counts.values()], labels=counts.keys())
plt.savefig('multilen_compare.png')


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
