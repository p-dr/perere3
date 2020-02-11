# description: Aligns RNASeq trimmed reads from SRA to S. mansoni genome database (ht2 formatted). Outputs SAM files.
# in: pardir/'trimmed_SRA_data' pardir/'genome_ht2'
# out: pardir/'alinhamentos/SRA_vs_genoma'

from glob import iglob
from pathlib import Path
from subprocess import call
from sys import argv
from utils import pardir, log
import multiprocessing as mp
from datetime import datetime

trimmed_data_dir = pardir/'trimmed_SRA_data'
out_dir = pardir/'alinhamentos/SRA_vs_genoma'
genome_prefix = str(pardir/'genome_ht2/sm_genome')

accs = set(Path(filepath).stem.split('_')[0] for filepath in iglob(str(trimmed_data_dir/'*')))


def align_acc(acc):
    input_prefix = str(trimmed_data_dir/acc)
    out_path = out_dir/acc

    if not out_path.exists() or '-r' in argv:
        print(f'Alinhando biblioteca {acc} com genoma.')
        t = datetime.now()
        call((f'hisat2 -x {genome_prefix} -1 {input_prefix}_1_paired* -2 {input_prefix}_2_paired* '
              f'--rna-strandness FR -S {str(out_path)+".sam"} -p {mp.cpu_count()}'), shell=True)
        dt = datetime.now() - t

        # sp.call(f'samtools sort -O BAM -no {str(out_path)+".bam"}')  # Will we use bam?

        print (f'Hisat2 terminou de rodar. Tempo decorrido: {dt}')
        log(f"Alinhamentos de {acc} com genoma concluídos e salvos em '{str(out_path)}'. Tempo decorrido: {dt}", __file__)        
    else:
        print (f"'{str(out_path)}' já existe. Pulando '{acc}'.")


def main():
    ## Deixar paralização por conta do Hisat2.
    # with mp.Pool() as pool:
    #     pool.map(align_acc, accs)
    log(f'Nova sessão de alinhamentos iniciada.')

    t = datetime.now()
    for acc in accs:
        align_acc(acc)
    dt = datetime.now() - t

    log(f'Sessão concluída em {dt}.')


if __name__ == '__main__':
    main()
