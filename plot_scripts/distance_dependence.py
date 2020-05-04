# in: d_path
import pandas as pd
# from statsmodels.nonparametric.smoothers_lowess import lowess
import utils as u
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append(str(u.scripts_dir))
from aggregate_data import outpath as d_path

d = pd.read_table(d_path)

plt.figure(figsize=(16, 9))
plt.yscale('log')
plt.xscale('log')
d = d[d.transcription != 0]
sns.regplot(d.distance, d.transcription, lowess=True, marker='.',
            scatter_kws={'alpha': .2})

cols = ('transcription', 'correlation', 'gene_transcription',
        'complement_transcription', 'complement_correlation', 'gene_complement_transcription',
       )

for col in cols:
    print(f'Plotting {col}...')
    plt.figure(figsize=(16, 9))
    plt.xlim(5e3, 5e4)
    plt.title(col)
    plt.yscale('log')
    plt.xscale('log')

    for label, subset in u.get_subsets(d).items():
        sns.regplot(subset.distance, subset[col],
                    lowess=True, label=label, scatter=False)
    plt.legend()

u.save_all_figs()
