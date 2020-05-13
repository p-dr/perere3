# description: Downloads RNASeq reads from SRA online database.
# in: (cloud)
# out: SRA_data_dir

import subprocess as sp
import multiprocessing as mp
from sys import argv
from glob import glob
import utils as u
from trim_SRA_data import trimmed_data_dir
from align_SRA_to_genome import out_dir as SRA_to_genome_alignments
import count_SRA_reads

n_cpu = mp.cpu_count()
SRA_data_dir = u.pardir/'SRA_data'
chunk_size = 2  # Maximum number of accession codes fetched in a run
# 4 as 4 files per acc, so 16 files for 16 cores. (2 because disk space is running out).
FINISHED_FLAG = 'finished_fetching'

# accs = ['ERR0228'+str(i) for i in range(72, 82)] + ['SRR922067', 'SRR922068']  # Aparently not stranded.
# accs = (['SRX74506'+str(i) for i in range(72, 79)] +
#         ['SRX74506'+str(i) for i in range(84, 87)])

# import accs from file
with open('seqs/SRAaccs.txt', 'r') as accs_file:
    accs = [acc.strip() for acc in accs_file.readlines()
            if not (acc.startswith('#') or acc == '\n')]


def fetch_acc(acc):
    existing_files = (
        *count_SRA_reads.out_dir.glob(acc+'*'),
        *SRA_to_genome_alignments.glob(acc+'*'),
        *trimmed_data_dir.glob(acc+'*'),
        *SRA_data_dir.glob(acc+'*'),
    )

    u.log(existing_files,'not existing', not bool(existing_files), 'redo', u.redo_flag)
    # if no outfile in outdir
    if not existing_files or u.redo_flag:
        u.log(f"Baixando '{acc}'...")
        try:
            sp.run(f'fasterq-dump --split-files -O {str(SRA_data_dir)} {acc} -e {n_cpu} -t /dev/shm -p', shell=True, check=True)
            # run(f'fastq-dump --gzip --split-files -I -O {str(SRA_data_dir)} {acc}', shell=True)

        except sp.CalledProcessError as err:
            u.log(err)
            raise err

    else:
        u.log(f"Pulando '{acc}' por existirem os seguintes arquivos:",
            *existing_files, sep='\n')
        return 'skipped'


def main():
    # You had better do this management in other place, thoug I understand this is the pipeline starting point.

    if u.args.clean:  # another flag should be created.
        leftover_files = (
            *SRA_to_genome_alignments.glob('*'),
            *trimmed_data_dir.glob('*'),
            *SRA_data_dir.glob('*'),
        )

        if leftover_files:
            u.log('Encountered following leftover files:',
                  *leftover_files, sep='\n\t', end='\n\n')

            while True:
                inp = input('Delete leftover files? [y/n]:')
                if inp in 'yY':
                    u.clean(*leftover_files)
                    break
                elif inp in 'nN':
                    raise RuntimeError('There are leftover files in the pipeline.')
                

    i = 0
    for acc in accs:
        if i == chunk_size: 
            u.log(f'Chunk downloaded ({chunk_size} files).')
            return
        if fetch_acc(acc) != 'skipped':
           i += 1 

    u.log('Finished downloading:', accs, sep=' ')
    return FINISHED_FLAG


if __name__ == '__main__':
    main()
