# description: Removes adapters and low quality ends from SRA raw data (trims it).
# in: pardir/'SRA_data'
# out: pardir/'trimmed_SRA_data'

from glob import iglob
from pathlib import Path
import subprocess as sp
from sys import argv
from utils import pardir, log, redo_flag
from wget import download
import multiprocessing as mp
from datetime import datetime


SRA_data_dir = pardir/'SRA_data'
trimmed_data_dir = pardir/'trimmed_SRA_data'
trimmed_data_dir.mkdir(exist_ok=True)

trimmomatic_source = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip'
trimmomatic_dir = pardir/'Trimmomatic-0.39'

if not trimmomatic_dir.exists():
    print('Trimmomatic não encontrado. Baixando...')
    download(trimmomatic_source, out=str(pardir))
    call(f'unzip {str(trimmomatic_dir)}.zip -d {str(pardir)}', shell=True)

trimmomatic_path = trimmomatic_dir/'trimmomatic-0.39.jar'
adapters_dir = trimmomatic_dir/'adapters'
adapter_path = str(adapters_dir/'TruSeq3-PE.fa')

EXT = '.fastq.gz'

accs = set(Path(filepath).stem.split('_')[0] for filepath in iglob(str(SRA_data_dir/'*')))

# Remember you chose _1 suffix for forward and _2 for reverse.
def trim_acc(acc):
    input_prefix = str(SRA_data_dir/acc)
    output_prefix = str(trimmed_data_dir/acc)
    output_sample_path = trimmed_data_dir/(acc+'_1_paired'+EXT)

    if not output_sample_path.exists() or redo_flag:
        t = datetime.now()
        try:
            output = sp.check_output((
                f'java -jar {str(trimmomatic_path)} PE -phred33 {input_prefix}_1{EXT} {input_prefix}_2{EXT} '
                f'{output_prefix}_1_paired{EXT} {output_prefix}_1_unpaired{EXT} '
                f'{output_prefix}_2_paired{EXT} {output_prefix}_2_unpaired{EXT} '
                f'ILLUMINACLIP:{adapter_path}:2:30:10 LEADING:3 TRAILING:3 SLIDIN'
                 'GWINDOW:4:15 MINLEN:36'), shell=True, universal_newlines=True)
            dt = t - datetime.now()

        except sp.CalledProcessError as err:
            log(f'ERROR: {err.returncode}\n OUT: {err.output}')

        else:
            log('STDOUT:\n' + output)

        log(f'{acc} trimado em {dt}.')

    else:
        print (f"'{str(output_sample_path)}' já existe. Pulando '{acc}'.")

with mp.Pool() as p:
    t = datetime.now()
    p.map(trim_acc, accs)
    dt = t-datetime.now()
    log(f'Sessão encerrada. Tempo total decorrido: {dt}')
