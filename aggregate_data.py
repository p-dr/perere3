# description: Plots scatter matrix comparing all new parameters and generate
# all_together_now.tsv compilation.
# in: pardir/'genome_annotation/head_genes_relations_unconsidering_sense.tsv'
# in: pardir/'genome_annotation/head_genes_correlations_unconsidering_sense.tsv'
# in: pardir/'genome_annotation/head_annotations.gff3'
# in: pardir/'counted_reads/aggregated_unconsidering_sense.tsv'
# in: pardir/'genome_annotation/gene_annotations.gff3'
# in: pardir/'genome_annotation/heads_repetitions.tsv'
# out: pardir/'genome_annotation/all_together_now.tsv'

import pandas as pd
import matplotlib.pyplot as plt
from utils import (pardir, save_all_figs, GFF3_COLUMNS,
                   unfold_gff, plot_box, show_flag)

outpath = pardir/'genome_annotation/all_together_now.tsv'

def reverse(s):
    return s[::-1]


# ##### RELATION FLAG WITH NEIGHBOR GENE
rel_data = pd.read_table(pardir/'genome_annotation' /
                         'head_genes_relations_unconsidering_sense.tsv')
rel_data.set_index('head_id', inplace=True)
rel_data.rename(columns={'gene_id': 'neighbor_gene',
                         'flag':'relative_position'}, inplace=True)

#############
used_genes = rel_data.neighbor_gene.unique()
#############

# ##### CORRELATION WITH NEIGHBOR GENE
corr_data = pd.read_table(pardir/'genome_annotation' /
                          'head_genes_correlations_unconsidering_sense.tsv')
corr_data.set_index('head_id', inplace=True)
corr_data.rename(columns={'gene_id': 'neighbor_gene'}, inplace=True)


# ##### HEAD DATA FROM ANNOTATION FILE
gff_data = pd.read_table(pardir/'genome_annotation/head_annotations.gff3',
                         header=None, names=GFF3_COLUMNS)
gff_data = unfold_gff(gff_data)
gff_data.set_index('gene_id', inplace=True)
gff_data.motherlength = gff_data.motherlength.astype(dtype='int')
gff_data.drop(['source', 'type', 'score', 'phase'], 1, inplace=True)

# ##### HEAD REPETITIONS (# OF INTERHEAD BLAST ALIGNMENTS)
heads_repetitions = pd.read_table(pardir/'genome_annotation/heads_repetitions.tsv', index_col='head_id')

# ##### READ COUNTS DATA
count_data = pd.read_table(pardir/'counted_reads' /
                           'aggregated_unconsidering_sense.tsv')
count_data = count_data.sum()[1:]
count_data.rename('transcription', inplace=True)
heads_count = count_data.loc[count_data.index.str.startswith('head')]
genes_count = count_data.loc[used_genes]

# ##### GENE ANNOTATIONS DATA
gene_gff = pd.read_table(pardir/'genome_annotation/gene_annotations.gff3',
                         header=None, names=GFF3_COLUMNS)
SELECTED_COLS = ['start', 'end', 'strand', 'length', 'Name']
gene_gff = unfold_gff(gene_gff)[SELECTED_COLS]
gene_gff = gene_gff.set_index('Name').loc[used_genes]
#gene_gff.drop(['source', 'type', 'score', 'phase', 'seqid'], 1, inplace=True)

# ##### COMPILE GENES DATA
genes_data = pd.concat([gene_gff, genes_count], 1)
genes_data.columns = ['gene_' + name for name in genes_data.columns]

# ##### ALL TOGETHER NOW!
data = pd.concat([gff_data, rel_data, corr_data.correlation, heads_count, heads_repetitions], 1, sort=False)
print(heads_repetitions)
print('='*100)
print(data)
print('='*100)
print(genes_data)
print('='*100)
data = data.merge(genes_data, left_on='neighbor_gene', right_index=True,
                  how='left')

# ##### FURTHER INFO APPENDS

# If there is no neighbor gene, 'same_strand' = False.
data['same_strand'] = (data.strand==data.gene_strand) & ~data.gene_strand.isna()

# ## Annotate stream.

# data['gene_stream'] = pd.np.nan
# data.loc[data.relative_position == 'olap', 'stream'] = 'olap'

# data.loc[((data.strand == '+') & (data.relative_position == 'gh')) | 
#          ((data.strand == '-') & (data.relative_position == 'hg')), 'gene_stream'] = 'gh'
# 
# data.loc[((data.strand == '+') & (data.relative_position == 'hg')) | 
#          ((data.strand == '-') & (data.relative_position == 'gh')), 'gene_stream'] = 'hg'

data['gene_stream'] = data.relative_position
mask = (data.gene_strand == '-') & data.relative_position.isin({'hg', 'gh'}) 
data.loc[mask, 'gene_stream'] = data.loc[mask, 'gene_stream'].apply(reverse)

data['head_stream'] = data.relative_position
mask = (data.strand == '-') & data.relative_position.isin({'hg', 'gh'}) 
data.loc[mask, 'head_stream'] = data.loc[mask, 'head_stream'].apply(reverse)


# ##### FINAL PRINTS

print('Quantidade de dados não faltantes:')
print((~data.isna()).sum())

print('\nDados com sondas relacionadas a genes mas sem correlação')
print((data.relative_position.isna() ^ data.correlation.isna()).sum())

print('''\nNotas: Só são exibidos os genes que foram relacionados a alguma head/sonda. Tem
      menos dados de genes porque tem algumas poucas sondas que não tinham
      nenhum gene no mesmo contig. Há bastante perda de dados quando se quer
      calcular a correlação, já que muitas sondas tiveram RPKM não-nulo em
      menos de 4 bibliotecas e foram então descartadas.''')


############ RESET INDEX AND SAVE ##################

data.index.name = 'head_id'
data = data.reset_index()

data.to_csv(outpath, sep='\t', index=False)
print('\nDados agregados salvos. Plotando...')

####################################################


# ###### PLOT ALL NUMERIC!
data.dropna(inplace=True)
data.drop('same_strand', 1, inplace=True)
data = data.infer_objects()
data = data.select_dtypes(exclude=['object'])

pd.plotting.scatter_matrix(data.drop(['end', 'start'], 1))
data.plot.scatter('start', 'transcription', alpha=.2)
data.plot.scatter('distance', 'transcription', alpha=.2)
plot_box(data[data.distance > 0], 'distance', 'transcription', bins=50)

genome_map = data.sort_values('start')[['start', 'transcription']]
genome_map['start'] = pd.cut(genome_map.start, 1e4)
genome_map = genome_map.groupby('start').sum()
le = len(genome_map.transcription)
plt.figure()
plt.pcolor(genome_map.transcription.values.reshape(int(le**.5), int(le**.5)))
plt.colorbar()
plt.figure()

plt.pcolor(data.corr(method='spearman'))
plt.figure()
data.transcription.hist()

if not show_flag:
    save_all_figs()
else:
    plt.show()
