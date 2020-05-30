# description: Aligns RNASeq trimmed reads from SRA to S. mansoni genome database (ht2 formatted). Outputs SAM files.
# in: trimmed_data_dir u.pardir/'genome_ht2'
# out: out_dir

from pathlib import Path
from sys import argv
from datetime import datetime
from count_SRA_reads import counted  # should not be here (will cause circular import after fixing import in count_SRA_reads)
import utils as u
import subprocess as sp
import multiprocessing as mp
from trim_SRA_data import trimmed_data_dir

genome_prefix = str(u.pardir/'genome_ht2/sm_genome')

out_dir = u.pardir/'alinhamentos/SRA_vs_genoma'
out_dir.mkdir(parents=True, exist_ok=True)


def align_acc(acc):
    input_prefix = str(trimmed_data_dir/acc)
    out_path = out_dir/(acc+".sam")

    if not counted(acc) or not out_path.exists() or '-r' in argv:
        try:
            print(f'Alinhando biblioteca {acc} com genoma.')
            t = datetime.now()
            hisat = sp.run((f'hisat2 -q -x {genome_prefix} -1 {input_prefix}_1_paired* -2 {input_prefix}_2_paired* '
                    f'--rna-strandness FR -S {str(out_path)} -p {mp.cpu_count()}'),
                    shell=True, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True, check=True)
            dt = datetime.now() - t

            u.log(hisat.stdout)
            # sp.call(f'samtools sort -O BAM -no {str(out_path)+".bam"}')  # Will we use bam?

            u.log(f"Alinhamentos de {acc} com genoma salvos em '{str(out_path)}'. Tempo decorrido: {dt}")        
            u.clean(*trimmed_data_dir.glob(acc+'*'))

        except sp.CalledProcessError as err:
            u.log(f'{acc} returned ERROR: {err.returncode} | OUT: {err.output}')
            #u.clean(out_path)
            raise err
            
    else:
        print (f"'{str(out_path)}' já existe ou ja foi contado. Pulando '{acc}'.")


def main():
    from trim_SRA_data import trimmed_data_dir
    accs = {p.stem.split('_')[0] for p in trimmed_data_dir.glob('*')}

    u.log(f'Nova sessão de alinhamentos iniciada.')
    u.log('Accession codes:', *accs, sep='\n\t')

    t = datetime.now()
    for acc in accs:
        align_acc(acc)
    dt = datetime.now() - t

    u.log(f'Sessão concluída em {dt}.')


if __name__ == '__main__':
    main()
