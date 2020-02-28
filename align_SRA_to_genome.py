# description: Aligns RNASeq trimmed reads from SRA to S. mansoni genome database (ht2 formatted). Outputs SAM files.
# in: pardir/'trimmed_SRA_data' pardir/'genome_ht2'
# out: pardir/'alinhamentos/SRA_vs_genoma'

from glob import iglob
from pathlib import Path
from sys import argv
from datetime import datetime
from utils import pardir, log, LOG_PATH
import subprocess as sp
import multiprocessing as mp

trimmed_data_dir = pardir/'trimmed_SRA_data'
genome_prefix = str(pardir/'genome_ht2/sm_genome')

out_dir = pardir/'alinhamentos/SRA_vs_genoma'
out_dir.mkdir(parents=True, exist_ok=True)

accs = set(Path(filepath).stem.split('_')[0] for filepath in iglob(str(trimmed_data_dir/'*')))


def align_acc(acc):
    input_prefix = str(trimmed_data_dir/acc)
    out_path = out_dir/acc

    if not out_path.exists() or '-r' in argv:
        print(f'Alinhando biblioteca {acc} com genoma.')
        t = datetime.now()
        sp.run((f'hisat2 -q -x {genome_prefix} -1 {input_prefix}_1_paired* -2 {input_prefix}_2_paired* '
                f'--rna-strandness FR -S {str(out_path)+".sam"} -p {mp.cpu_count()}'),
                shell=True, stdout=LOG_PATH.open('a'), stderr=sp.STDOUT)
        dt = datetime.now() - t

        # sp.call(f'samtools sort -O BAM -no {str(out_path)+".bam"}')  # Will we use bam?

        print (f'Hisat2 terminou de rodar. Tempo decorrido: {dt}')
        log(f"Alinhamentos de {acc} com genoma concluídos e salvos em '{str(out_path)}'. Tempo decorrido: {dt}")        
    else:
        print (f"'{str(out_path)}' já existe. Pulando '{acc}'.")


def main():
    log(f'Nova sessão de alinhamentos iniciada.')

    t = datetime.now()
    for acc in accs:
        align_acc(acc)
    dt = datetime.now() - t

    log(f'Sessão concluída em {dt}.')


if __name__ == '__main__':
    main()
