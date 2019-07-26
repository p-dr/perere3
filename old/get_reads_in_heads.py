from utils import read_tsv, pardir, redo_flag, overlaps, verbose, log, progress_flag
from re import findall
from time import time
# NÃO ESTÁ SE LEVANDO EM CONSIDERAÇÃO OS SENTIDOS DAS FITAS.

samdir = pardir/'alinhamentos'/'SRA_vs_genoma'
outdir = pardir/'alinhamentos'/'SRA_vs_genoma_on_heads'

    
def len_from_cigar(cigar):
    return sum([int(i) for i in findall(r'\d+', cigar)])


if __name__ == '__main__':
    
    heads = read_tsv(pardir/'genome_annotation'/'head_annotations.gff3', header=None)
    heads['id'] = heads[8].apply(lambda s: s.strip('gene_id=').split(';')[0])
    heads = heads.iloc[:, [0, 3, 4, -4, -1]]
    heads.columns = ['chrom', 'start', 'end', 'sense', 'id']
    head_groups = heads.groupby('chrom')
    
    for sam_path in samdir.glob('*.sam'):

        outpath = outdir/(sam_path.stem + '.tsv')

        if outpath.exists() and not redo_flag:
            print(f"'{str(outpath)}' existe. Use -r se quiser sobrescrever.")
            continue

        print(f"\nTrabalhando com arquivo '{str(sam_path)}'.")
        #================== get n lines of sam file ===================#
        if progress_flag:
            print('Calculando número total de linhas...')
            with sam_path.open('rb') as sam_file:
                total_len = sum(1 for i in sam_file)

        print(f"Escrevendo arquivo '{str(outpath)}'...\n")
        outfile = outpath.open('w')

        outfile.write('\t'.join(['chrom', 'head_id',
                                 'head_start', 'head_end',
                                 'read_start', 'read_cigar'])+'\n')

        sam_df = read_tsv(sam_path, comment='@', header=None, names=range(22), chunksize=1)

        count = 0
        t0 = time()
        for sam_row in sam_df:
            sam_row = sam_row[[2, 3, 5]]
            sam_row.columns = ['chrom', 'start', 'cigar']

            if progress_flag:
                print(f'\rLinha SAM: {count} | Progresso: {count/total_len*100:.4f}% | Tempo restante: {(time()-t0)*(total_len/count-1):.4f}s', end='')
            count += 1

            sam_row = sam_row.iloc[0, :]

            if sam_row.chrom == '*':
                continue
            
            try:
                head_group = head_groups.get_group(sam_row.chrom)

            except KeyError:
                if verbose:
                    print(f'\nNão há heads no cromossomo {sam_row.chrom}.')
                continue
            
            for _, head_row in head_group.iterrows():

                if overlaps((sam_row.start, sam_row.start + len_from_cigar(sam_row.cigar)),
                            (head_row.start, head_row.end)):

                    outfile.write('\t'.join([head_row.chrom,
                                             head_row.id, str(head_row.start),
                                             str(head_row.end), str(sam_row.start),
                                             sam_row.cigar])+'\n')
            
        outfile.close()
        log(f"O SAM '{sam_path.stem}' acabou de ser filtrado.", __file__)
        print('Concluído.')
