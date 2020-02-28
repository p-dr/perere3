# description:  Generates filtered_perere3_vs_genoma, which is perere3_vs_genoma without alignments also and better found for SR3. i.e. Removes each perere alignment overlapped with SR3 alignment of best score. Also writes filtered out line's indices.
# in: pardir/'alinhamentos/perere3_vs_genome.bl'
# in: pardir/'alinhamentos/sr3_vs_genome.bl'
# out: pardir/'alinhamentos/perere3_vs_genome_sr3_filtered.bl'
# out: pardir/'alinhamentos/perere3_vs_genome_discarded.tsv'

import pandas as pd
import numpy as np
from utils import overlaps, pardir, verbose, BL_COLUMNS, prinf, safe_open
import multiprocessing as mp
from tqdm import tqdm
from datetime import datetime

# Usar biopython para blastear?
filtered_outpath = pardir/'alinhamentos/perere3_vs_genome_sr3_filtered.bl'
discarded_outpath = pardir/'alinhamentos/perere3_vs_genome_discarded.tsv'
filtered_outfile = safe_open(filtered_outpath, exist_ok=False)
discarded_outfile = safe_open(discarded_outpath, exist_ok=False)

perere3_inpath = pardir/'alinhamentos/perere3_vs_genome.bl'
sr3_inpath = pardir/'alinhamentos/sr3_vs_genome.bl'

n_cpu = mp.cpu_count()

#================== LER E FILTRAR ALINHAMENTOS ==================#

print('Lendo resultados do Blast...', end=' ')
perere3_vs_genoma = pd.read_table(perere3_inpath,
                             header=None, names=BL_COLUMNS)
sr3_vs_genoma = pd.read_table(sr3_inpath,
                         header=None, names=BL_COLUMNS)
print('Resultados lidos.')

# Sort positions
# for data in (perere3_vs_genoma, sr3_vs_genoma):
#     data.sort_values('sstart', inplace=True)
#     data.reset_index(drop=True, inplace=True)


def cartesian_product(left, right):
    return left.assign(tmp=0).merge(right.assign(tmp=0), on='tmp').drop('tmp', 1)


def parse_product(progress_pos, product):
    chunk_discarded = pd.DataFrame()
    tqdm.pandas(position=progress_pos+2, leave=False, desc=f'Thread {progress_pos:2}')
    overlapping_mask = product[['sstart_x', 'send_x', 'sstart_y', 'send_y']].progress_apply(overlaps, axis=1, raw=True)

    if not overlapping_mask.empty:
        overlapping = product.loc[overlapping_mask]
        chunk_discarded = chunk_discarded.append(overlapping)
        return chunk_discarded


def main():
 
    print('Filtrando alinhamentos em que o SR3 é melhor...')
    discarded = pd.DataFrame()
    filtered_perere3_vs_genoma = perere3_vs_genoma.copy()

    p_groups = perere3_vs_genoma.groupby('saccver')
    s_groups = sr3_vs_genoma.groupby('saccver') ## agrupar muda index???????????

    print('Iterando para cada scaffold no genoma e para cada perere3 no scaffold.')
    for p_group_name in tqdm(p_groups.groups, desc='Scaffolds'):
        s_group = s_groups.get_group(p_group_name)
        p_group = p_groups.get_group(p_group_name)

        prinf('Combinando DataFrames...', end='\r')
        product = cartesian_product(p_group[['sstart', 'send', 'bitscore']].reset_index(),
                                    s_group[['sstart', 'send', 'bitscore']].reset_index())

        # discard when perere3 aligns better
        prinf('Filtrando por bitscore do SR3...', end='\r')
        product = product.loc[product.bitscore_x < product.bitscore_y,
                              ]

        if product.empty:
            continue

        prinf('Subdividindo produto...         ', end='\r')
        product_chunks = np.array_split(product, n_cpu)

        prinf('Procurando sobreposições...', end='\r')
        with mp.Pool() as pool:
            chunks_discarded = pool.starmap(parse_product, enumerate(product_chunks))

        group_discarded = pd.concat(chunks_discarded)
        discarded = discarded.append(group_discarded)
        # print(discarded[~discarded['index'].isin(filtered_perere3_vs_genoma.index)])
        # print(filtered_perere3_vs_genoma.loc[discarded['index'].unique()])


    filtered_perere3_vs_genoma.drop(discarded['index'].values, inplace=True)

    print(f'\nFiltragem concluída. {len(discarded)} alinhamentos removidos.')

    print(f"Escrevendo alinhamentos filtrados do perere3 em '{str(filtered_outpath)}'...", end=' ')
    filtered_perere3_vs_genoma.to_csv(filtered_outfile, sep='\t', index=False)
    print('Arquivo escrito.')

    print(f"Escrevendo posições das linhas removidas de '{str(perere3_inpath)}' em '{str(discarded_outpath)}'...", end=' ')

    discarded.columns = pd.MultiIndex.from_product([('perere3', 'sr3'), ('index', 'sstart', 'ssend', 'bitscore')])
    discarded.to_csv(discarded_outfile, sep='\t', index=False)
    print('Arquivo escrito.')

    return filtered_perere3_vs_genoma, discarded


if __name__ == '__main__':
    filtered_perere3_vs_genoma, discarded = main()
