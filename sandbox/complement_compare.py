import sys
sys.path.append('..')

import utils as u
import pandas as pd
import matplotlib.pyplot as plt

counts_dir = u.pardir/'scripts/test/counts'

print('WARNING: Be aware we are imposing .csv extension for count outputs, which should be .tsv in the future.')
forward_reverse_pairs = [(p, p.parent/(p.stem + '_complement.csv'))
                        for p in counts_dir.glob('*')
                        if 'complement' not in p.stem]

# for pair in forward_reverse_pairs:
#     print(pair)
# exit()

n_pairs = len(forward_reverse_pairs)
plt.figure(figsize=(4*n_pairs, 8))

for i, pair in enumerate(forward_reverse_pairs):
    genes_path, complements_path = pair
    title = genes_path.stem
    u.print_header(title)

    plt.subplot(1, n_pairs, i+1)
    plt.title(title)
    genes, complements = [pd.read_table(p, header=None, skipfooter=5, engine='python')
                          for p in (genes_path, complements_path)]

    print(pd.concat([genes.describe(), complements.describe()], 1))
    print('\nNon-zeros:', (genes[1] != 0).sum(), (complements[1] != 0).sum())
    u.box_compare(genes[1], complements[1], labels=['genes', 'complements'])

plt.tight_layout()
plt.savefig('graficos/complement_compare.png')
