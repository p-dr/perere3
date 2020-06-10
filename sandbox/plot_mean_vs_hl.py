import pandas as pd
import matplotlib.pyplot as plt

for f in ('multiheadlen_transcription.tsv', 'multiheadlen_complement_transcription.tsv'):
    d = pd.read_table(f)
    d.plot.scatter('headlen', 'mean', figsize=(16,9))
    plt.savefig('graficos/' + f + '.png')
