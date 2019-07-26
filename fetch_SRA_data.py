# description: Downloads RNASeq reads from SRA online database.
# in: (cloud)
# out: pardir/'SRA_data'

from subprocess import run
from sys import argv
from glob import glob
from utils import redo_flag as redo, pardir

SRA_data_dir = pardir/'SRA_data'

accs = ['ERR0228'+str(i) for i in range(72, 82)] + ['SRR922067', 'SRR922068']
accs = accs[int(argv[1]):int(argv[2])]

    
for acc in accs:
    output_wc = str(SRA_data_dir/(acc+'*'))
    
    if not glob(output_wc) or redo:
        print(f"Baixando '{acc}'...")
        run(f'fastq-dump --gzip --split-files -I -O {str(SRA_data_dir)} {acc}', shell=True)
        print ('Comando de download executado.')

    else:
        print (f"'{output_wc}' já existe. Pulando '{acc}'.")
