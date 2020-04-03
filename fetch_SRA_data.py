# description: Downloads RNASeq reads from SRA online database.
# in: (cloud)
# out: pardir/'SRA_data'

from subprocess import run
import multiprocessing as mp
from sys import argv
from glob import glob
from utils import redo_flag, pardir, log
from trim_SRA_data import trimmed_data_dir
import count_SRA_reads

n_cpu = mp.cpu_count()
SRA_data_dir = pardir/'SRA_data'

# accs = ['ERR0228'+str(i) for i in range(72, 82)] + ['SRR922067', 'SRR922068']  # Aparently not stranded.
accs = (['SRX74506'+str(i) for i in range(72, 79)] +
        ['SRX74506'+str(i) for i in range(84, 87)])

# accs = [
#     'SRX7450686',
#     'SRX7450685',
#     'SRX7450672',
# ]


def fetch_acc(acc):
    existing_files = (
        *count_SRA_reads.out_dir.glob(acc+'*'),
        *trimmed_data_dir.glob(acc+'*'),
        *SRA_data_dir.glob(acc+'*'),
    )

    # if no outfile in outdir
    if not existing_files or redo_flag:
        log(f"Baixando '{acc}'...")
        exit()
        run(f'fasterq-dump --split-files -O {str(SRA_data_dir)} {acc} -e {n_cpu} -t /dev/shm -p', shell=True)
        # run(f'fastq-dump --gzip --split-files -I -O {str(SRA_data_dir)} {acc}', shell=True)
        log('Comando de download executado.')

    else:
        log(f"Pulando '{acc}' por existirem os seguintes arquivos:",
            *existing_files, sep='\n')


def main():
    for acc in accs:
        fetch_acc(acc)


if __name__ == '__main__':
    main()
