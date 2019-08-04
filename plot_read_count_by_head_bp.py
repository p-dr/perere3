# description: Plots how many read were aligned at each head's bp.
# in: pardir/'reads_by_heads_bp'
# out: 
# plot: 

from utils import pardir, redo_flag, verbose, argv
from pickle import load
import matplotlib.pyplot as plt

indir = pardir/'reads_by_heads_bp'

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

    #normalizar
    max_count = max(counts.values())
    if 'norm' in argv:
        complete_counts = [i/max_count for i in complete_counts]

    if max_count > THRESHOLD:
        plt.subplot(fig_rows, 1, contador%fig_rows+1)
        plt.plot(complete_counts, label=head)
        if contador%fig_rows == fig_rows-1:
            plt.show()
        contador += 1
        
plt.title('Contagem de reads por base das heads')
plt.xlabel('Posição na head')
plt.ylabel('Contagem')
if contador <= 10:
    plt.legend()
plt.show()
print('Pronto.')
