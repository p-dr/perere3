# description: Calculates expression (read count) correlation between heads and its nearest gene(s). When only "reation" is said, I mean head -> nearest gene(s) mapping.
# in: pardir/'genome_annotation/head_genes_relations.tsv' pardir/'counted_reads' pardir/'genome_annotation/gene_annotations.gff3' pardir/'genome_annotation/head_annotations.gff3'
# out: pardir/'genome_annotation/head_genes_correlations.tsv' pardir/'counted_reads/aggregated.tsv'

from utils import read_tsv, pardir, show_flag, redo_flag, parse_gff_attributes,
prinf, GFF3_COLUMNS, save_all_figs
import pandas as pd
from matplotlib import pyplot as plt
from sys import argv

if '--sem_sentido' in argv:
    nosense_flag = '_unconsidering_sense'

else:
    nosense_flag = ''


outpath = pardir/f'genome_annotation/head_genes_correlations{nosense_flag}.tsv'
out_aggregated_counts = pardir/f'counted_reads/aggregated{nosense_flag}.tsv'

if outpath.exists() and not (redo_flag or show_flag):
    print(f"'{str(outpath)}' já existe, nada será feito.")
    exit()

if redo_flag:
    outfile = (outpath).open('w')

print('Buscando comprimentos de genes e heads...')
gene_attibutes = read_tsv(pardir/'genome_annotation/gene_annotations.gff3', names=GFF3_COLUMNS, usecols=['attributes'])['attributes']
head_attibutes = read_tsv(pardir/'genome_annotation/head_annotations.gff3', names=GFF3_COLUMNS, usecols=['attributes'])['attributes']
gene_lengths = parse_gff_attributes(gene_attibutes)['length']
head_lengths = parse_gff_attributes(head_attibutes)['length']

# df.infer_objects() can't do '3' -> 3.  (-_-)
lengths = pd.concat([head_lengths, gene_lengths]).astype(int)

print('Concluído. Lendo arquivo de relações...')
relations = read_tsv(pardir/f'genome_annotation/head_genes_relations{nosense_flag}.tsv')

if redo_flag:
    outfile.write('\t'.join(relations.columns) + '\tcorrelation\n')

print('Conclúido. Agregando contagens em um DataFrame...')

# Read counts for each feature
counts = pd.DataFrame()

# for each SRA library
for count_path in (pardir/'counted_reads').glob('*.csv'):

    lib_name = count_path.stem.split('_')[0]
    
    count = read_tsv(count_path, header=None,
                     names=['feature', lib_name],
                     index_col='feature')

    # This block does not run if count is being generated by counter script
    if not count.empty:
        
        ### NORMALIZE FOR TOTAL LIBRARY SIZE (total_count):
        # total_count = counts sum
        # + __no_feature
        # + __ambiguous
        # + __too_low_aQual
        # (without) __not_aligned
        # + __alignment_not_unique
        total_count = count.iloc[:-2].sum() + count.iloc[-1]
        count /= total_count/1e6
        
        ### NORMALIZE FOR FEATURE LENGTH (IN KBP):
        count[lib_name] /= lengths/1000
        count.dropna(inplace=True)

        # count final meaning:
        # read count by library size (in million reads) by feature length (in bp)
        # we can call that something like  feature's "expression coeficient" in that library.
        counts = counts.combine_first(count)
        

counts = counts.T
counts.index.name = 'biblioteca'
#counts.reset_index(inplace=True)
if redo_flag:
    counts.to_csv(out_aggregated_counts, sep='\t')

print('Concluído. Calculando correlações...')


gridx, gridy = 5, 8
plt.rcParams['font.size'] = 4
plt.rcParams['figure.figsize'] = (16, 9)
i = 0

for _, relation_row in relations.iterrows():

    hid = relation_row.head_id
    gid = relation_row.gene_id
    
    if hid not in counts:
        prinf(f'WARNING: {hid} não presente nas contagens, só nas relações head-gene. Talvez a contagem deva ser refeita.')
        continue
    head_col = counts[hid]
    gene_col = counts[gid]

    n_non_zero = (head_col.astype(bool) & gene_col.astype(bool)).sum()
    # if number of significative points (head and gene counts != 0) is less than 4
    if n_non_zero < 4:
        prinf('\nCorrelação descartada:\n', counts[[gid, hid]])
        continue

    prinf('\nCorrelação plotada:\n', pd.concat([counts[[gid, hid]], head_col.astype(bool) & gene_col.astype(bool)], 1))

    corr = head_col.corr(gene_col)

    if show_flag:
        gridi = i % (gridx*gridy) + 1
        plt.subplot(gridx, gridy, gridi)
        plt.plot(head_col, gene_col, 'o', label=corr)
        plt.legend()
        plt.title('_'.join(list(relation_row.astype(str)) + [str(n_non_zero)]))

        if gridi == gridx*gridy:
            plt.tight_layout()
            save_all_figs()
            plt.show()

        i += 1
        
    elif redo_flag:
        outfile.write('\t'.join([*relation_row.astype(str), str(corr)])+'\n')

if redo_flag:
    outfile.close()