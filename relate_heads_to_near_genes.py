# description: Maps each head to its nearest (or overlapped) gene.
# in: pardir/'genome_annotation/head_annotations.gff3' pardir/'genome_annotation/gene_annotations.gff3'
# out: pardir/'genome_annotation/head_genes_relations.tsv'
# pardir/'genome_annotation/head_genes_relations_unconsidering_sense.tsv'

from utils import pardir, read_tsv, overlaps, redo_flag, parse_gff_attributes, prinf, GFF3_COLUMNS
from sys import argv

# Achar genes sobrepostos (ou o mais próximo se não se sobrepuserem), às heads pra correlacionar as expressões
# em um script posterior. Gera tsv em outpath.
# LEGENDA: gene está à dir(direita)/esq(esquerda)/olap(overlappado com) a head.

GFF_COLS_SUBSET = ['seqid', 'start', 'end', 'strand', 'attributes']

# ESTAMOS CONSIDERANDO SENTIDO POR DEFAULT (--sem_sentido para não considerar)
# sequid é o nome do cromossomo (contig)
if '--sem_sentido' in argv:
    COLS_TO_GROUP = 'seqid'
    nosense_flag = '_unconsidering_sense'

else:
    COLS_TO_GROUP = ['seqid', 'strand']
    nosense_flag = ''

outpath = pardir/f'genome_annotation/head_genes_relations{nosense_flag}.tsv'

# Segurança para não sobrescrever como eu fiz agora >.<
if outpath.exists() and not redo_flag:
    print(f"Arquivo '{str(outpath)}' já existe, nada será feito. Use '-r' se quiser sobrescrever.")
    exit()

outfile = outpath.open('w')

heads = read_tsv(pardir/'genome_annotation/head_annotations.gff3', names=GFF3_COLUMNS, usecols=GFF_COLS_SUBSET)
genes = read_tsv(pardir/'genome_annotation/gene_annotations.gff3', names=GFF3_COLUMNS, usecols=GFF_COLS_SUBSET)
heads['id'] = parse_gff_attributes(heads.attributes).index
genes['id'] = parse_gff_attributes(genes.attributes).index
head_groups = heads.groupby(COLS_TO_GROUP)

# ### alterar um pouquinho as duplicatas muahahah
genes.loc[genes.duplicated('start', 'last'), 'start'] += 1
genes.loc[genes.duplicated('end', 'last'), 'end'] += 1

gene_groups = genes.groupby(COLS_TO_GROUP)


if __name__ == '__main__':
    #write header
    outfile.write('\t'.join(['head_id', 'gene_id', 'flag', 'distance'])+'\n')

    for head_group_name, head_group in head_groups:
        try:
            gene_group = gene_groups.get_group(head_group_name)

        except KeyError:
            prinf('Não há nenhum gene no cromossomo. As heads abaixo são "desgenadas".')
            prinf(head_group)
            continue

        for _, head_row in head_group.iterrows():

            for _, gene_row in gene_group.iterrows():
                if overlaps((gene_row.start, gene_row.end),
                            (head_row.start, head_row.end)):
                    flag = 'olap'
                    chosen_gene_id = gene_row.id
                    distance = 0
                    break

            # if none overlaps
            else:
                # b f de forward ou backfill
                distance = 9e99
                for meth in ['b', 'f']:
                    fill_back = meth == 'b'
                    gene_pos = ['end', 'start'][fill_back]
                    head_pos = ['end', 'start'][not fill_back]
                    indexed_group = gene_group.sort_values(gene_pos).set_index(gene_pos)

                    try:
                        near_gene_index = indexed_group.index.get_loc(head_row[head_pos], meth+'fill')
                        gene_row = gene_group.iloc[near_gene_index]
                        new_distance = abs(gene_row[gene_pos] - head_row[head_pos])

                        if new_distance < distance:
                            chosen_gene_id = gene_row.id
                            flag = ["esq", "dir"][fill_back]
                            distance = new_distance

                    except KeyError:
                        prinf(f'Não há gene à {flag} de {head_row.id}.')

            outfile.write('\t'.join([head_row.id, chosen_gene_id, flag, str(distance)])+'\n')

    print('Concluído.')

outfile.close()
