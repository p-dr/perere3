# description: Removes adapters and low quality ends from SRA raw data (trims it).
# in: SRA_data_dir
# out: trimmed_data_dir

import subprocess as sp
from sys import argv
import utils as u
from wget import download
import multiprocessing as mp
from datetime import datetime
import count_SRA_reads
from align_SRA_to_genome import out_dir as SRA_to_genome_alignments


SRA_data_dir = u.pardir/'SRA_data'
trimmed_data_dir = u.pardir/'trimmed_SRA_data'

trimmomatic_source = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip'
trimmomatic_dir = u.pardir/'Trimmomatic-0.39'

trimmomatic_path = trimmomatic_dir/'trimmomatic-0.39.jar'
adapters_dir = trimmomatic_dir/'adapters'
adapter_path = str(adapters_dir/'TruSeq3-PE.fa')


# Remember you chose _1 suffix for forward and _2 for reverse.
def trim_acc(acc):
    inp1, inp2 = sorted(list(SRA_data_dir.glob(acc + '*')))
    ext = inp2.suffix
    output_prefix = str(trimmed_data_dir/acc)

    existing_files = (
        *count_SRA_reads.out_dir.glob(acc+'*'),
        *SRA_to_genome_alignments.glob(acc+'*'),
        *trimmed_data_dir.glob(acc+'*'),
    )

    # if outfiles are not in outdir
    if not existing_files or u.redo_flag:
        t = datetime.now()
        try:
            command = (
                f'java -jar {str(trimmomatic_path)} PE -threads {mp.cpu_count()} -phred33 {str(inp1)} {str(inp2)} '
                f'{output_prefix}_1_paired{ext} {output_prefix}_1_unpaired{ext} '
                f'{output_prefix}_2_paired{ext} {output_prefix}_2_unpaired{ext} '
                f'ILLUMINACLIP:{adapter_path}:2:30:10 LEADING:3 TRAILING:3 SLIDIN'
                 'GWINDOW:4:15 MINLEN:36')
            output = sp.check_output(command, shell=True, universal_newlines=True)
            dt = datetime.now() - t
            u.log(f'{acc} trimado em {dt}.')
            u.clean(inp1, inp2)

        except sp.CalledProcessError as err:
            u.log(f'ERROR: {err.returncode}\n OUT: {err.output}')
            #u.clean(*trimmed_data_dir.glob(acc+'*'))
            raise err


    else:
        u.log(f"Pulando '{acc}' por existirem os seguintes arquivos:",
            *existing_files, sep='\n')


def main():
    trimmed_data_dir.mkdir(exist_ok=True)

    if not trimmomatic_dir.exists():
        print('Trimmomatic não encontrado. Baixando...')
        download(trimmomatic_source, out=str(u.pardir))
        sp.run(f'unzip {str(trimmomatic_dir)}.zip -d {str(u.pardir)}', shell=True)
    accs = {filepath.stem.split('_')[0] for filepath in SRA_data_dir.glob('*')}

    t = datetime.now()
    list(map(trim_acc, accs))
    dt = datetime.now() - t
    u.log(f'Sessão encerrada. Tempo total decorrido: {dt}')


if __name__ == '__main__':
    main()
