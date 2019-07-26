# description: Maps each head to its nearest (or overlapped) gene.
# in: pardir/'genome_annotation/head_annotations.gff3' pardir/'genome_annotation/gene_annotations.gff3'
# out: pardir/'genome_annotation/head_genes_relations.tsv'

from utils import pardir, read_tsv, overlaps, verbose, redo_flag

# Achar genes sobrepostos (ou o mais próximo se não se sobrepuserem), às heads pra correlacionar as expressões
# em um script posterior. Gera tsv em outpath.
# LEGENDA: gene está à dir(direita)/esq(esquerda)/olap(overlappado com) a head.

outpath = pardir/'genome_annotation/head_genes_relations.tsv'
outfile = outpath.open('w')

# Segurança para não sobrescrever como eu fiz agora >.<
if outpath.exists() and not redo_flag:
    print(f"Arquivo '{str(outpath)}' já existe, nada será feito. Use '-r' se quiser sobrescrever.")
    exit()

    
GFF_COLS = ['chrom', 'start', 'end', 'sense', 'id']

# ESTAMOS CONSIDERANDO SENTIDO
heads = read_tsv(pardir/'genome_annotation/head_annotations.gff3', header=None)
heads.drop([1, 2, 5, 7], inplace=True, axis=1)
heads[8] = heads[8].apply(lambda s: s.strip('gene_id=').split(';')[0])
heads.columns = GFF_COLS
head_groups = heads.groupby(['chrom', 'sense'])

genes = read_tsv(pardir/'genome_annotation/gene_annotations.gff3', header=None)
genes.drop([1, 2, 5, 7], inplace=True, axis=1)
genes[8] = genes[8].apply(lambda s: s.strip('gene_id=').split(';')[0])
genes.columns = GFF_COLS

# ### alterar um pouquinho as duplicatas muahahah
genes.loc[genes.duplicated('start', 'last'), 'start'] += 1
genes.loc[genes.duplicated('end', 'last'), 'end'] += 1

gene_groups = genes.groupby(['chrom', 'sense'])


if __name__ == '__main__':
    #write header
    outfile.write('\t'.join(['head_id', 'gene_id', 'flag', 'distance'])+'\n')
    
    for head_group_name, head_group in head_groups:
        try:
            gene_group = gene_groups.get_group(head_group_name)

        except KeyError:
            if verbose:
                print('Não há nenhum gene no cromossomo. As heads abaixo são "desgenadas".')
                print(head_group)
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
                        if verbose:
                            print(f'Não há gene à {flag} de {head_row.id}.')

            outfile.write('\t'.join([head_row.id, chosen_gene_id, flag, distance])+'\n')

    print('Concluído.')

outfile.close()
