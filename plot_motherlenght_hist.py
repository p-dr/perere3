# description: Plots heads' motherlength distribution.
# in: pardir/'genome_annotation/heads_motherlength.tsv'
# out: 
# plot: 

from pandas import read_csv
from matplotlib import pyplot
from utils import pardir, save_all_figs

motherlengths = read_csv(str(pardir/'genome_annotation' /
                             'heads_motherlength.tsv'),
                         sep='\t',
                         header=None,
                         names=['head_id', 'motherlength'])
motherlengths.hist(bins=50, figsize=(9, 4.8))
pyplot.xlabel('Comprimento-mãe (bp)')
pyplot.ylabel('Frequência')
pyplot.title('Distribuição de comprimentos-mãe das cópias de Pererê-3')
save_all_figs()

pyplot.show()
