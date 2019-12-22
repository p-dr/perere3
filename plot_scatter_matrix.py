# description: Plots scatter matrix comparing all new parameters and generate
# all_together_now.tsv compilation.
# in: pardir/'genome_annotation/head_genes_relations_unconsidering_sense.tsv'
# in: pardir/'genome_annotation/head_genes_correlations_unconsidering_sense.tsv'
# in: pardir/'genome_annotation/head_annotations.gff3'
# in: pardir/'counted_reads/aggregated_unconsidering_sense.tsv'
# in: pardir/'genome_annotation/gene_annotations.gff3'
# out: pardir/'genome_annotation/all_together_now.tsv'

import pandas as pd
import matplotlib.pyplot as plt
from utils import (pardir, save_all_figs, GFF3_COLUMNS,
                   unfold_gff, plot_box, show_flag)

outpath = pardir/'genome_annotation/all_together_now.tsv'

# ##### RELATION FLAG WITH NEIGHBOR GENE
rel_data = pd.read_table(pardir/'genome_annotation' /
                         'head_genes_relations_unconsidering_sense.tsv')
rel_data.set_index('head_id', inplace=True)
rel_data.rename(columns={'gene_id': 'neighbor_gene'}, inplace=True)

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
gene_gff = unfold_gff(gene_gff)
gene_gff = gene_gff.set_index('gene_id').loc[used_genes]
gene_gff.drop(['source', 'type', 'score', 'phase', 'seqid'], 1, inplace=True)

# ##### COMPILE GENES DATA
genes_data = pd.concat([gene_gff, genes_count], 1)
genes_data.columns = ['gene_' + name for name in genes_data.columns]

# ##### ALL TOGETHER NOW!
data = pd.concat([gff_data, rel_data, corr_data.correlation, heads_count], 1)
data = data.merge(genes_data, left_on='neighbor_gene', right_index=True,
                  how='left')

print('Dados não faltantes:')
print((~data.isna()).sum())

data['same_strand'] = data.strand == data.gene_strand

print('\nDados com sondas relacionadas a genes mas sem correlação')
print((data.flag.isna() ^ data.correlation.isna()).sum())

print('''\nNotas: Só são exibidos os genes que foram relacionados a alguma head/sonda. Tem
      menos dados de genes porque tem algumas poucas sondas que não tinham
      nenhum gene no mesmo contig. Há bastante perda de dados quando se quer
      calcular a correlação, já que muitas sondas tiveram RPKM não-nulo em
      menos de 4 bibliotecas e foram então descartadas.''')

##############################
data.to_csv(outpath, sep='\t')
print('\nDados agregados salvos.')
##############################

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
plt.show()
