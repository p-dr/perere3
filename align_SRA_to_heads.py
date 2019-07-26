# description: Aligns RNASeq trimmed reads from SRA to heads database (ht2 formatted). Outputs SAM files.
# in: pardir/'trimmed_SRA_data' pardir/'heads_ht2/heads'
# out: pardir/'alinhamentos/SRA_vs_heads'

from glob import iglob
from pathlib import Path
from subprocess import call
from utils import pardir, log, redo_flag

# It was written before as 'trimmed_with_perere_SRA_data'.
trimmed_data_dir = pardir/'trimmed_SRA_data'
out_dir = pardir/'alinhamentos/SRA_vs_heads'

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
        call(f'hisat2 -x {heads_prefix} -1 {input_prefix}_1_paired* -2 {input_prefix}_2_paired* -S {str(out_path)}', shell=True)
        print ('Hisat2 terminou de rodar.')
        log(f"Alinhamentos de {acc} com heads concluídos e salvos em '{str(out_path)}'.", __file__)
        
    else:
        print (f"'{str(out_path)}' já existe. Pulando '{acc}'.")

print ('Pronto.')
