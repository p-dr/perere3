import sys
sys.path.append('..')

import utils as u
import pandas as pd
import matplotlib.pyplot as plt

counts_dir = u.pardir/'scripts/test/counts'

genes, complements = [pd.read_table(p, header=None, skipfooter=5)
                      for p in counts_dir.glob('*')]

print('GENES', genes.head(1))
print('COMPLEMENTS', complements.head(1))

print(pd.concat([genes.describe(), complements.describe()], 1))

plt.figure(figsize=(4, 8))
u.multibox_compare([genes[1], complements[1]], labels=['genes', 'complements'])
plt.tight_layout()
plt.savefig('graficos/complement_compare.png')
plt.show()

