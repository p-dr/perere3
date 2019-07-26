# description: Calculates expression (read count) correlation between heads and its nearest gene(s). When only "reation" is said, I mean head -> nearest gene(s) mapping.
# in: pardir/'genome_annotation/head_genes_relations.tsv' pardir/'counted_reads'
# out: pardir/'genome_annotation/head_genes_correlations.tsv'

from utils import read_tsv, pardir, show_flag, redo_flag, verbose
import pandas as pd
from matplotlib import pyplot as plt

outpath = pardir/'genome_annotation/head_genes_correlations.tsv'

if outpath.exists() and not (redo_flag or show_flag):
    print(f"'{str(outpath)}' já existe, nada será feito.")
    exit()

if redo_flag:
    outfile = (outpath).open('w')

print('Lendo arquivo de relações...')
relations = read_tsv(pardir/'genome_annotation/head_genes_relations.tsv', header=None, names=['head_id', 'gene_id', 'flag'])

print('Conclúido. Agregando contagens em um DataFrame...')

counts = pd.DataFrame()

for count_path in (pardir/'counted_reads').glob('*.csv'):
    count = read_tsv(count_path, header=None,
                     names=['feature', count_path.stem.split('_')[0]],
                     index_col='feature')

    # if you are running this script while counting reads
    if not count.empty:
        total_count = count.iloc[:-2].sum() + count.iloc[-1]
        counts = counts.combine_first(count/total_count*1e6)

counts = counts.T
counts.index.name = 'biblioteca'
counts.reset_index(inplace=True)
print('Concluído. Calculando correlações...')


gridx, gridy = 5, 8
plt.rcParams['font.size'] = 4
plt.rcParams['figure.figsize'] = (16, 9)
i = 0

for _, relation_row in relations.iterrows():

    head_col = counts[relation_row.head_id]
    gene_col = counts[relation_row.gene_id]

    if head_col.astype(bool).sum() < 4:
        if verbose:
            print('\nHead descartada:\n', head_col)
        continue

    corr = head_col.corr(gene_col)

    if show_flag:
        gridi = i % (gridx*gridy) + 1
        plt.subplot(gridx, gridy, gridi)
        plt.plot(head_col, gene_col, 'o', label=corr)
        plt.legend()
        plt.title('_'.join(tuple(relation_row)))

        if gridi == gridx*gridy:
            plt.tight_layout()
            plt.show()

        i += 1
    else:
        outfile.write('\t'.join([*relation_row, str(corr)])+'\n')

outfile.close()
