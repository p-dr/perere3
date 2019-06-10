#!/bin/python

# Generates PARENT_DIR/alinhamentos/heads_vs_heads.bl if '-r' provided
# and plots N number of alignments versus number of sequences aligned N
# times to graficos/heads_vs_heads.png and graficos/heads_vs_heads_log.png

from PARENT_DIR import PARENT_DIR
from os.path import exists
from pandas import read_csv
from sys import argv
from matplotlib import pyplot as plt

COLUMNS = 'qaccver saccver pident mismatch gapopen qstart qend sstart send evalue bitscore'


#======================== CRIAR HEADS_VS_HEADS ========================#

heads_path = f'{PARENT_DIR}/seqs/heads.fa'
out_path = f'{PARENT_DIR}/alinhamentos/heads_vs_heads.bl'
out_not_exists = not exists(out_path)

if out_not_exists:
    print(f'{out_path} não existe.')

if out_not_exists or '-r' in argv:
    print('Procurando alinhamentos entre heads...')

    # Remember you are using megablast.
    run(f"blastn -query {heads_path} -subject {heads_path} -outfmt '6 {COLUMNS}' -out {out_path} -evalue 1e-10", shell=True)
    print(f'Alinhamentos salvos em {out_path}.\n')


#======================== LER ARQUIVO ========================#

print('Lendo alinhamentos heads_vs_heads...')
heads_vs_heads = read_csv(out_path, sep='\\s+', header=None, names=COLUMNS.split())
print(f"'{out_path}' lido.")


#======================== HISTOGRAMA ========================#

hist = (heads_vs_heads['qaccver'].value_counts()-1).value_counts()

for log in (False, True):

    print('Construindo histograma...')

    plt.rcParams['axes.titlesize'] = 24
    plt.rcParams['axes.labelsize'] = 20

    hist.plot('bar', grid=False, logy=log)

    plt.xlabel('Número de repetições')
    plt.ylabel('Número de sequências')
    plt.title('Número de sequências "head" que se repetem N vezes graficado contra o pŕoprio N'.upper())

    hist_path = f'{PARENT_DIR}/graficos/heads_vs_heads{["","_log"][log]}.png'
    plt.savefig(hist_path)
    print(f"Histograma salvo em '{hist_path}'.")
    
    if '-s' in argv:
        plt.show()
        plt.close()
