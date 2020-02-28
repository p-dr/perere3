# description: Maps each head to its nearest (or overlapped) gene.
# in: pardir/'genome_annotation/head_annotations.gff3'
# in: pardir/'genome_annotation/gene_annotations.gff3'
# out: pardir/'genome_annotation/head_genes_relations.tsv'

from utils import (pardir, overlaps,
                   redo_flag, parse_gff_attributes,
                   prinf, GFF3_COLUMNS, safe_open)
import pandas as pd
import numpy as np
from sys import argv
from tqdm import tqdm
import multiprocessing as mp

# Achar genes sobrepostos (ou o mais próximo se não se sobrepuserem), às heads pra correlacionar as expressões
# em um script posterior. Gera tsv em outpath.
# LEGENDA: 'gh' ou 'hg' representam a posição relativa entre cópia e NG ('hg'
# == cópia antes do gene vizinho), considerado a fita + e o sentido de escrita.

GFF_COLS_SUBSET = ['seqid', 'start', 'end', 'strand', 'attributes']

# sequid é o nome do cromossomo (contig)
COLS_TO_GROUP = 'seqid'
outpath = pardir/f'genome_annotation/head_genes_relations.tsv'
outfile = safe_open(outpath, exist_ok='exit')
use_multiprocessing = '--multiprocess' in argv
n_cpu = mp.cpu_count()

heads = pd.read_table(pardir/'genome_annotation/head_annotations.gff3', names=GFF3_COLUMNS, usecols=GFF_COLS_SUBSET)
genes = pd.read_table(pardir/'genome_annotation/gene_annotations.gff3', names=GFF3_COLUMNS, usecols=GFF_COLS_SUBSET)
heads['id'] = parse_gff_attributes(heads.attributes).index
genes['id'] = parse_gff_attributes(genes.attributes).index
head_groups = heads.groupby(COLS_TO_GROUP)
gene_groups = genes.groupby(COLS_TO_GROUP)


def parse_head_row(head_row, gene_group):
    global outfile

    parse_head_row.last_args = head_row, gene_group, outfile

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
        distance = 9e99  # create Inf() class for that?
        for meth in ['b', 'f']:
            fill_back = meth == 'b'
            gene_pos = ['end', 'start'][fill_back]
            head_pos = ['end', 'start'][not fill_back]

            if gene_group[gene_pos].duplicated().sum():
                raise pd.core.indexes.base.InvalidIndexError(
                    'You appear to have duplicated gene positions. '
                    'Indexes won\'t work properly.')

            indexed_group = gene_group.sort_values(gene_pos).set_index(gene_pos)

            try:
                near_gene_index = indexed_group.index.get_loc(head_row[head_pos], meth+'fill')
                gene_row = gene_group.iloc[near_gene_index]
                new_distance = abs(gene_row[gene_pos] - head_row[head_pos])

                if new_distance < distance:
                    chosen_gene_id = gene_row.id
                    flag = ["gh", "hg"][fill_back]
                    distance = new_distance

            except KeyError as e:
                prinf(f'\nNão há gene {meth.upper()} de {head_row.id}. (B=atrás, F=à frente)\n')

    print('\t'.join([head_row.id, chosen_gene_id, flag, str(distance)])+'\n')
    outfile.write('\t'.join([head_row.id, chosen_gene_id, flag, str(distance)])+'\n')


def parse_chunk(df, gene_group, n):
    global n_cpu
    tqdm.pandas(position=n+2, desc=str(n), leave=False)
    return df.progress_apply(parse_head_row, axis=1, args=[gene_group])


def main():
    # write header
    outfile.write('\t'.join(['head_id', 'gene_id', 'flag', 'distance'])+'\n')
    print('Iterate for each contig and for each head in contig.')

    for head_group_name, head_group in tqdm(head_groups, desc='Contigs'):
        try:
            gene_group = gene_groups.get_group(head_group_name)
            chunks = np.array_split(head_group, n_cpu)

            with mp.Pool() as pool:
                pool.starmap(parse_chunk, ((c, gene_group, cn) for cn, c in enumerate(chunks)))

        except KeyError:
            prinf('Não há nenhum gene no cromossomo. As heads abaixo são "desgenadas".')
            prinf(head_group)
            continue


    print(f'\nConcluído. Relações salvas em {str(outpath)}.')


if __name__ == '__main__':
    main()

outfile.close()
