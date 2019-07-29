# description: (CLARIFICATION NEEDED) Relates Perere3complete-to-genome alignments to heads. Adds lenght parameter to outfile.
# in: pardir/'alinhamentos/filtered_perere3complete_vs_genoma.bl' pardir/'genome_annotation/head_annotations.gff3'
# out: pardir/'genome_annotation/head_alignment_length.tsv' pardir/'alinhamentos/perere3complete_vs_geno-\nma_antisenso_invertida.bl'
# flags: unknown_purpose

from pandas import read_csv
from utils import pardir, verbose
from align_seqs_to_genome import COLUMNS
from gen_heads import HEAD_LEN

OUTPATH = pardir/'genome_annotation/head_alignment_length.tsv'


aligned_perere = read_csv(str(pardir/'alinhamentos/filtered_perere3complete_vs_genoma.bl'),
                          sep='\t', header=None, names=COLUMNS.split()).sort_values('send')  # Esse sort não devia ser feito após lidar com o sentido da fita (abaixo)?
head_annotations = read_csv(str(pardir/'genome_annotation/head_annotations.gff3'),
                            sep='\t', header=None, names=['acc', 'origin', 'kind',
                                                          'start', 'end', 'dot1',
                                                          'sense', 'dot2', 'attr'])


#========== LIDAR COM SENTIDO DA FITA ===========#
# Reconhecer linhas do arquivo BLAST que correspondem a alinhamentos na fita
# complementar, sinalizado pelo BLAST com ssend < sstart e põe-se os índices
# sstart e ssend na ordem crescente, visto que sstart será usado posteriorm-
# ente para ordenamento (foi?).

OUT_INVERTIDO = pardir/'alinhamentos/perere3complete_vs_genoma_antisenso_invertida.bl'
print('Lidando com sentido da fita...')
per_len = len(aligned_perere.index)
count = 0
for i, row in aligned_perere.iterrows():
    if row['sstart'] > row['send']:
        aligned_perere.loc[i, 'send'], aligned_perere.loc[i, 'sstart'] = row['sstart'], row['send']
        print(f'Sentido negativo encontrado. Progresso: {count/per_len*100:.4f}%', end='\r')
    count += 1
print('Valores trocados.')

print(aligned_perere)
#================================================#


head_annotations.set_index('start', inplace=True)
aligned_perere.set_index('send', inplace=True)

### Agrupar por cromossomo
# annotation são as head_annotations
aligned_groups = aligned_perere.groupby(['saccver'])
annotation_groups = head_annotations.groupby(['acc'])

al_len = len(aligned_perere)
an_len = len(head_annotations)

# annotations é bem menor. (define qual loop é o de dentro)(?)
# print('al', al_len, 'an', an_len)

SEARCH_LEN = 5
head_len_map = {}
count = 0

for acc in annotation_groups.groups:
    for start, annotation in annotation_groups.get_group(acc).iterrows():
        target_slice = slice(start-SEARCH_LEN, annotation['end']-HEAD_LEN+SEARCH_LEN)
        #target_slice = slice(0, 1)
        match = aligned_groups.get_group(acc).loc[target_slice]
        if not match.empty:
            if len(match.index) > 1:
                print(f'WARNING: Ambiguidade encontrada na região {target_slice}: \n{match}')

            length = match['length'].values[0]

            if verbose:
                print(start, annotation['end'], f"{annotation['attr']} -> lenght={length}")

            # GRAVAR MATCH!
            head_len_map[int(annotation['attr'].split('=')[-1].strip('head'))] = length
            
            print(f'Há matches sendo encontrados. Progresso: {count/an_len*100:.4f}%', end='\r')
        count+=1

#print(head_len_map.items(), sep='\n')

print(f"Escrevendo '{str(OUTPATH)}'...")
with OUTPATH.open('w') as outfile:
    for headid, l in sorted(head_len_map.items()):
        outfile.write(f'{headid}\t{l}\n')

print('Arquivo escrito.')
