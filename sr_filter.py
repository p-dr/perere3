# description:  Generates filtered_perere3_vs_genoma, which is perere3_vs_genoma without alignments also and better found for SR3. i.e. Removes each perere alignment overlapped with SR3 alignment of best score. Also writes filtered out line's indices.
# in: pardir/'alinhamentos/perere3_vs_genome.bl'
# in: pardir/'alinhamentos/sr3_vs_genome.bl'
# out: pardir/'alinhamentos/filtered_perere3_vs_genome.bl'
# out: pardir/'alinhamentos/perere3_indices_filtrados.csv'

from pandas import read_csv
from utils import overlaps, pardir, verbose, BL_COLUMNS
from tqdm import tqdm

# Usar biopython para blastear?


#================== LER E FILTRAR ALINHAMENTOS ==================#

print('Lendo resultados do Blast...', end=' ')
perere3_vs_genoma = read_csv(pardir/'alinhamentos/perere3_vs_genome.bl',
                             header=None, names=BL_COLUMNS, sep='\\s+')
sr3_vs_genoma = read_csv(pardir/'alinhamentos/sr3_vs_genome.bl',
                         header=None, names=BL_COLUMNS, sep='\\s+')
print('Resultados lidos.')

perere3_vs_genoma = perere3_vs_genoma.sort_values(['saccver', 'sstart'])
sr3_vs_genoma = sr3_vs_genoma.sort_values(['saccver', 'sstart'])

sr3_grouped = sr3_vs_genoma.groupby('saccver')
perere3_grouped = perere3_vs_genoma.groupby('saccver')

for group_name in perere3_grouped.groups:
    p_group = perere3_grouped.get_group(group_name)
    s_group = sr3_grouped.get_group(group_name)
    print(p_group, s_group)

print('\nFiltragem concluída')

filtered_outpath = pardir/'alinhamentos/filtered_perere3complete_vs_genoma.bl'
print(f"Escrevendo alinhamentos filtrados do perere3 em '{str(filtered_outpath)}'...")
filtered_perere3_vs_genoma.to_csv(str(filtered_outpath), sep='\t', index=False)
print('Arquivo escrito.')

discarded_outpath = pardir/'alinhamentos/perere3complete_indices_filtrados.csv'
print(f"Escrevendo posições das linhas removidas de 'perere3complete_vs_genoma.bl' em '{str(discarded_outpath)}'...")
discarded_outpath.open('w').write('\n'.join(discarded))
print('Arquivo escrito.')
