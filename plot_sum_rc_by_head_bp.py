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

heads_data = unfold_gff(heads_data)
print(heads_data)

print('Agregando dicionários...')
total_count_dic = {}

for count_dic_path in indir.glob('*'):

    with count_dic_path.open('rb') as count_dic_file:
        count_dic = load(count_dic_file)

    for head in count_dic.keys():
        head_total_count = sum(count_dic[head].values())
        for pos, count in count_dic[head].items():
            total_count_dic[pos] = total_count_dic.get(pos, 0) + count/head_total_count

print('Feito. Completando posições faltantes e plotando...')

complete_counts = [total_count_dic.get(pos, 0) for pos in
                   range(max(total_count_dic.keys())+1)]
print('Posições foram completadas. Plotando...')

motherlength = motherlengths.loc[head, 'motherlength']

axesy = plt.axis()[2:]
n_heads = len(heads_data)
mask_count = [0]*1000
count = 0

### FIND AND PLOT REPETIVIVE REGIONS
def plot_repetitive():
    try:
        for irow, row in heads_data.iterrows():
            seqid, hstart, hend = row[['seqid', 'start', 'end']]

            for j, mask in mask_map.get_group(seqid).iterrows():
                # if mask does not overlaps
                if mask.end < hstart or mask.start > hend:
                    continue
                else:
                    olap_mask = (max(mask.start - hstart, 0), min(mask.end - hstart, 1000))
                    for i in range(*olap_mask):
                        mask_count[i] += 1
                    # plt.fill_between(olap_mask, *axesy, alpha=.1)
                    count += 1

            print(irow, '/', n_heads, 'plotted:', count)

    except KeyboardInterrupt:
        pass

    return mask_count


# NORMALIZE
complete_counts = [c/max(complete_counts) for c in complete_counts]
#mask_count = [c/max(mask_count) for c in mask_count]

plt.plot(complete_counts, label='Read count')
#plt.plot(mask_count, label='Masked count')
plt.title('Contagem de reads por base das heads')
plt.xlabel('Posição na head')
plt.ylabel('Contagem')
save_all_figs()
plt.show()
print('Pronto.')
