from subprocess import run
from sys import argv
from pathlib import Path
from glob import glob

pardir = Path(__file__).resolve().parents[1]
SRA_data_dir = pardir/'SRA_data'

accs = ['ERR0228'+str(i) for i in range(72, 82)] + ['SRR922067', 'SRR922068']

redo = '-r' in argv
if redo:
    argv.remove('-r')

if len(argv) > 1:
    accs = accs[int(argv[1]):int(argv[2])]

    
for acc in accs:
    output_wc = str(SRA_data_dir/(acc+'*'))
    
    if not glob(output_wc) or redo:
        print(f"Baixando '{acc}'...")
        run(f'fastq-dump --gzip --split-files -I -O {str(SRA_data_dir)} {acc}', shell=True)
        print ('Comando de download executado.')

    else:
        print (f"'{output_wc}' já existe. Pulando '{acc}'.")
