# description: Aligns (BLAST) heads between themselves. Plots the number of head sequences aligned N times to other heads as a function of N itself.
# in: pardir/'seqs/heads.fa'
# out: pardir/'alinhamentos/heads_vs_heads.bl'
# out: pardir/'alinhamentos/heads_repetitions.tsv'
# plot: 

from utils import pardir, safe_open, show_flag, save_all_figs
from subprocess import run
from os.path import exists
from pandas import read_csv
from sys import argv
from matplotlib import pyplot as plt

COLUMNS = 'qaccver saccver pident mismatch gapopen qstart qend sstart send length evalue bitscore'


#======================== CRIAR HEADS_VS_HEADS ========================#

heads_path = pardir/'seqs/heads.fa'
out_path = pardir/'alinhamentos/heads_vs_heads.bl'
out_file = safe_open(out_path)  # Check out_path.
repetitions_outpath = pardir/'genome_annotation/heads_repetitions.tsv'

if out_file is not None:
    # Remember you are using megablast.
    run(f"blastn -task 'megablast' -query '{heads_path}' -subject '{heads_path}' -outfmt '6 {COLUMNS}' -out '{out_path}' -evalue 1e-10", shell=True)
    print(f'Alinhamentos salvos em {out_path}.\n')


#======================== LER ARQUIVO ========================#

print('Lendo alinhamentos heads_vs_heads...')
heads_vs_heads = read_csv(out_path, sep='\\s+', header=None, names=COLUMNS.split())
print(f"'{out_path}' lido.")


#======================== HISTOGRAMA ========================#

counts = (heads_vs_heads['qaccver'].value_counts()-1)
counts_df = counts.to_frame().reset_index()
counts_df.columns = ['head_id', 'repetitions']
counts_df.to_csv(repetitions_outpath, sep='\t', index=False)

hist = counts.value_counts()
print(hist)
soma = hist.sum()
print('soma:', soma, 'fração de 0 repetições:',
      hist[0]/soma, 'fração com até 1 rep.', hist[0:2].sum()/soma)

for log in (False, True):

    print(f'Construindo histograma, log={log}')
    plt.figure(figsize=(11, 4.8))

    # ax = counts.hist(bins=200, figsize=(4, 4.8))
    plt.bar(hist.index, hist, log=log)

    plt.xlabel('Número de repetições')
    plt.ylabel('Número de sequências')
    plt.title('Quantidade de cópias repetidas')

    print("Histograma salvo.")

    if show_flag:
        plt.show()

if show_flag:
    save_all_figs()

plt.close()
