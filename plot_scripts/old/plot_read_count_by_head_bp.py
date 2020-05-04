# description: Plots how many read were aligned at each head's bp.
# in: pardir/'reads_by_heads_bp' pardir/'genome_annotation/heads_motherlength.tsv'
# out: 
# plot: 

from utils import pardir, redo_flag, verbose, argv, overlaps, unfold_gff, GFF3_COLUMNS, save_all_figs
from pickle import load
from scipy.signal import find_peaks
import pandas as pd
import matplotlib.pyplot as plt

# Minimum motherlength value to be plotted.
ML_THRESH = 3000

mask_map = pd.read_csv(pardir/'genome_annotation/mask_map.tsv',
                       sep='\t').groupby('seqid')

indir = pardir/'reads_by_heads_bp'
motherlengths = pd.read_csv(pardir/'genome_annotation/heads_motherlength.tsv',
                            sep='\t', header=None, names=['head_id', 'motherlength'],
                            index_col='head_id')

heads_data = pd.read_csv(pardir/'genome_annotation/head_annotations.gff3',
                         names=GFF3_COLUMNS, sep='\t', usecols=['seqid', 'start', 'end',
                                                      'strand', 'attributes'])

heads_data = unfold_gff(heads_data).set_index('gene_id')
print(heads_data)

print('Agregando dicionários...')
total_count_dic = {}

for count_dic_path in indir.glob('*'):

    with count_dic_path.open('rb') as count_dic_file:
        count_dic = load(count_dic_file)

    for head in count_dic.keys():
        total_count_dic[head] = total_count_dic.get(head, {})
        for pos, count in count_dic[head].items():
            total_count_dic[head][pos] = total_count_dic[head].get(pos, 0) + count

print('Feito. Completando posições faltantes e plotando...')
contador = 0

try:
    THRESHOLD = int(argv[-1])
except ValueError:
    THRESHOLD = 0

fig_rows = 5

for head, counts in total_count_dic.items():
    complete_counts = [counts.get(pos, 0) for pos in range(max(counts.keys())+1)]

    # Normalizar
    max_count = max(counts.values())
    if 'norm' in argv:
        complete_counts = [i/max_count for i in complete_counts]

    motherlength = motherlengths.loc[head, 'motherlength']

    if max_count > THRESHOLD and motherlength > ML_THRESH:

        plt.subplot(fig_rows, 1, contador%fig_rows+1)
        plt.legend()

        plt.plot(complete_counts, label=f'{head} (ml={motherlength})')
        axesy = plt.axis()[2:]

        ### FIND AND PLOT REPETIVIVE REGIONS
        seqid, hstart, hend = heads_data.loc[head, ['seqid', 'start', 'end']]

        masks = []
        for j, mask in mask_map.get_group(seqid).iterrows():
            # if mask does not overlaps
            if mask.end < hstart or mask.start > hend:
                continue
            else:
                olap_mask = (mask.start - hstart, mask.end - hstart)
                print(*zip(olap_mask, axesy))
                plt.fill_between(olap_mask, *axesy, alpha=.5)


        if contador%fig_rows == fig_rows-1:
            save_all_figs()
            plt.show()
        contador += 1

plt.title('Contagem de reads por base das heads')
plt.xlabel('Posição na head')
plt.ylabel('Contagem')
if contador <= 10:
    plt.legend()
save_all_figs()
plt.show()
print('Pronto.')
