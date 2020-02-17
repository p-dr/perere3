# description: Aligns RNASeq trimmed reads from SRA to heads database (ht2 formatted). Outputs SAM files.
# in: pardir/'trimmed_SRA_data' pardir/'heads_ht2'
# out: pardir/'alinhamentos/SRA_vs_heads'

from glob import iglob
from pathlib import Path
from subprocess import run
from utils import pardir, log, redo_flag, LOG_PATH
from multiprocessing import cpu_count

# Future: align from cloud with --sra-acc option; only one command: one can input multiple files to Hisat2.
# It was written before as 'trimmed_with_perere_SRA_data'.
trimmed_data_dir = pardir/'trimmed_SRA_data'
out_dir = pardir/'alinhamentos/SRA_vs_heads'
log_file = LOG_PATH.open('a')
ncpu = cpu_count()

# to use with --un-conc
#unpaired_out_dir = out_dir/'despareados'

heads_prefix = str(pardir/'heads_ht2/heads')

accs = set(Path(filepath).stem.split('_')[0] for filepath in iglob(str(trimmed_data_dir/'*')))

for acc in accs:
    input_prefix = str(trimmed_data_dir/acc)
    out_path = out_dir/(acc+'.sam')
    #unpaired_out_path = unpaired_out_dir/(acc+'.sam')
    
    
    if not out_path.exists() or redo_flag:
        print (f"Alinhando {acc} com heads...")
        run((f'hisat2 -q -x {heads_prefix} -1 {input_prefix}_1_paired* '
             f'-2 {input_prefix}_2_paired* -S {str(out_path)} --no-unal --rna-strandness FR --threads {ncpu}'),
             shell=True, stout=log_file, stderr=log_file)
        print ('Hisat2 terminou de rodar.')
        log(f"Alinhamentos de {acc} com heads concluídos e salvos em '{str(out_path)}'.")
        
    else:
        print (f"'{str(out_path)}' já existe. Pulando '{acc}'.")

print ('Pronto.')
log_file.close()
