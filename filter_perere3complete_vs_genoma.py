#!/bin/python

# Generates filtered_perere3_vs_genoma, which is perere3_vs_genoma
# without alignments also found for SR3, i.e. overlapped
# with those of the latter.

from pandas import read_csv
from utils import overlaps, pardir, verbose
from time import time, gmtime, strftime
from align_seqs_to_genome import COLUMNS

# Usar biopython para blastear?


#================== LER E FILTRAR ALINHAMENTOS ==================#

print('Lendo resultados do Blast...')
perere3_vs_genoma = read_csv(pardir/f'alinhamentos/perere3complete_vs_genoma.bl', header=None, names=COLUMNS.split(), sep='\\s+')
sr3_vs_genoma = read_csv(pardir/f'alinhamentos/sr3complete_vs_genoma.bl', header=None, names=COLUMNS.split(), sep='\\s+')
print('Resultados lidos.')

for data in (perere3_vs_genoma, sr3_vs_genoma):
    data.sort_values('sstart', inplace=True)
    data.reset_index(drop=True, inplace=True)

print('Filtrando alinhamentos em que o SR3 é melhor...\n')
total_len = len(perere3_vs_genoma)
filtered_perere3_vs_genoma = perere3_vs_genoma.copy()

sr3_scaffold_groups = sr3_vs_genoma.groupby(['saccver'])
start_time = time()
count = 0
discarded = []
pos = 0


for i, p in perere3_vs_genoma.iterrows():

    count += 1

    if p['saccver'] in sr3_scaffold_groups.groups:
        for j, s in sr3_scaffold_groups.get_group(p['saccver']).loc[pos:].iterrows():


            if overlaps((p['sstart'], p['send']), (s['sstart'], s['send'])):

                # print('Alinhamentos coincidentes encontrados:')
                # print(f"{p['sstart']}-{p['send']}\n{s['sstart']}-{s['send']}")

                if s['bitscore'] > p['bitscore']:

                    # print('SR3 tem mais bitscore:')
                    # print(f"SR3 = {s['bitscore']}\nPer = {p['bitscore']}")
                    pos = j

                    filtered_perere3_vs_genoma.drop(i, inplace=True)
                    # print(f'Alinhamento {i} descartado.\n')
                    discarded.append(str(i))

                    if verbose:
                        elapsed = time() - start_time
                        elapsed_str = strftime('%H:%M:%S', gmtime(elapsed))
                        remaining_str = strftime('%H:%M:%S', gmtime(elapsed*total_len/count))
                        print(f'{elapsed_str}/{remaining_str} ({count/total_len*100:.4f}%) {pos}', end='\r')

                    break

                # else:
                #     print('Perere3 tem mais bitscore:')
                #     print(f"SR3 = {s['bitscore']}\nPer = {p['bitscore']}")
                #     print(f'Alinhamento {i} mantido.\n')


print('\nFiltragem concluída')

filtered_outpath = pardir/'alinhamentos/filtered_perere3complete_vs_genoma.bl'
print(f"Escrevendo alinhamentos filtrados do perere3 em '{str(filtered_outpath)}'...")
filtered_perere3_vs_genoma.to_csv(str(filtered_outpath), sep='\t', index=False)
print('Arquivo escrito.')

discarded_outpath = pardir/'scripts/perere3complete_indices_filtrados.csv'
print(f"Escrevendo posições das linhas removidas de 'perere3complete_vs_genoma.bl' em '{str(discarded_outpath)}'...")
discarded_outpath.open('w').write('\n'.join(discarded))
print('Arquivo escrito.')
