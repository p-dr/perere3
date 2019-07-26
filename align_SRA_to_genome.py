# description: Aligns RNASeq trimmed reads from SRA to S. mansoni genome database (ht2 formatted). Outputs SAM files.
# in: pardir/'trimmed_SRA_data' pardir/'genome_ht2/smgenome'
# out: pardir/'alinhamentos/SRA_vs_genoma'

from glob import iglob
from pathlib import Path
from subprocess import call
from sys import argv
from utils import pardir

trimmed_data_dir = pardir/'trimmed_SRA_data'
out_dir = pardir/'alinhamentos/SRA_vs_genoma'
genome_prefix = str(pardir/'genome_ht2/smgenome')

accs = set(Path(filepath).stem.split('_')[0] for filepath in iglob(str(trimmed_data_dir/'*')))

for acc in accs:
    input_prefix = str(trimmed_data_dir/acc)
    out_path = out_dir/(acc+'.sam')

    if not out_path.exists() or '-r' in argv:
        print (f"Alinhando {acc} com genoma...")
        call(f'hisat2 -x {genome_prefix} -1 {input_prefix}_1_paired* -2 {input_prefix}_2_paired* -S {str(out_path)}', shell=True)
        print ('Hisat2 terminou de rodar.')
        log(f"Alinhamentos de {acc} com genoma concluídos e salvos em '{str(out_path)}'.", __file__)        
    else:
        print (f"'{str(out_path)}' já existe. Pulando '{acc}'.")

print ('Pronto.')
