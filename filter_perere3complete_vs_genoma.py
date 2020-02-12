# description:  Generates filtered_perere3_vs_genome, which is perere3_vs_genome without alignments also and better found for SR3. i.e. Removes each perere alignment overlapped with SR3 alignment of best score. Also writes filtered out line's indices.
# in: pardir/'alinhamentos/perere3_vs_genome.bl' pardir/'alinhamentos/sr3_vs_genome.bl'
# out: pardir/'alinhamentos/filtered_perere3_vs_genome.bl' pardir/'alinhamentos/perere3_indices_filtrados.csv'

from pandas import read_csv
from utils import overlaps, pardir, verbose, BL_COLUMNS
from tqdm import tqdm

# Usar biopython para blastear?

tol = 4000
filtered_outpath = pardir/'alinhamentos/filtered_perere3_vs_genome_simple.bl'
discarded_outpath = pardir/'alinhamentos/perere3_indices_filtrados_simple.csv'

#================== LER E FILTRAR ALINHAMENTOS ==================#

print('Lendo resultados do Blast...', end=' ')
perere3_vs_genome = read_csv(pardir/'alinhamentos/perere3_vs_genome.bl', header=None, names=BL_COLUMNS, sep='\\s+')
sr3_vs_genome = read_csv(pardir/'alinhamentos/sr3_vs_genome.bl', header=None, names=BL_COLUMNS, sep='\\s+')
print('Resultados lidos.')

filtered_perere3_vs_genome = perere3_vs_genome.copy()

for data in (perere3_vs_genome, sr3_vs_genome):
    data.sort_values(['saccver', 'sstart'], inplace=True)

p_groups = perere3_vs_genome.groupby('saccver')
s_groups = sr3_vs_genome.groupby('saccver')

# for group_name in p_groups.groups:
#     p_group = p_group.get_group(group_name)
#     s_group = s_group.get_group(group_name)
#     pos = 0 
# 
#     for pi, p_row in p_group.iterrows():
#         for si, s_row in s_group.iloc[pos:].iterrows():
#             if overlaps(p_row.sstart, p_row.send,
#                         s_row.sstart, s_row.send):
#                  pos = si
                 
print('Filtrando alinhamentos em que o SR3 é melhor...\n')

for group_name in tqdm(p_groups.groups):

    p_group = p_groups.get_group(group_name)
    s_group = s_groups.get_group(group_name).set_index('sstart')
    print(s_group)

    for pi, p_row in tqdm(list(p_group.iterrows())):
        s_index = s_group.index.get_loc(p_row.sstart, method='nearest')
        s_row = s_group.loc[s_index]
        if overlaps(p_row.sstart, p_row.send,
                    s_row.sstart, s_row.send):

            if s_row.bitscore >= p_row.bitscore:
                filtered_perere3_vs_genome.drop(p_row.index, inplace=True)


print('\nFiltragem concluída')

print(f"Escrevendo alinhamentos filtrados do perere3 em '{str(filtered_outpath)}'...")
filtered_perere3_vs_genome.to_csv(str(filtered_outpath), sep='\t', index=False)
print('Arquivo escrito.')

print(f"Escrevendo posições das linhas removidas de 'perere3_vs_genome.bl' em '{str(discarded_outpath)}'...")
discarded_outpath.open('w').write('\n'.join(discarded))
print('Arquivo escrito.')
