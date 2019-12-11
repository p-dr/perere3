# description: Plots heads' motherlength distribution.
# in: pardir/'genome_annotation/heads_motherlength.tsv'
# out: 
# plot: 

from pandas import read_csv
from matplotlib import pyplot
from utils import pardir, save_all_figs, show_flag

motherlengths = read_csv(str(pardir/'genome_annotation' /
                             'heads_motherlength.tsv'),
                         sep='\t',
                         header=None,
                         names=['head_id', 'motherlength'])
motherlengths.hist(bins=50, figsize=(9, 4.8))
pyplot.xlabel('Comprimento (bp)')
pyplot.ylabel('Frequência')
pyplot.title('Distribuição de comprimentos das regiões conservadas das cópias de Perere-3')
save_all_figs()

if show_flag:
    pyplot.show()
