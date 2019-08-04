# description: Removes adapters and low quality ends from SRA raw data (trims it).
# in: pardir/'SRA_data' pardir/'trimming_adapters/TruSeq3-PE.fa'
# out: pardir/'trimmed_SRA_data'

from glob import iglob
from pathlib import Path
from subprocess import call
from sys import argv
from utils import pardir, log


SRA_data_dir = pardir/'SRA_data'
trimmed_data_dir = pardir/'trimmed_SRA_data'

adapters_dir = pardir/'trimming_adapters'
adapter_path = str(adapters_dir/'TruSeq3-PE.fa')

EXT = '.fastq.gz'

accs = set(Path(filepath).stem.split('_')[0] for filepath in iglob(str(SRA_data_dir/'*')))

# Remember you chose _1 suffix for forward and _2 for reverse.
for acc in accs:
    input_prefix = str(SRA_data_dir/acc)
    output_prefix = str(trimmed_data_dir/acc)
    output_sample_path = trimmed_data_dir/(acc+'_1_paired'+EXT)

    if not output_sample_path.exists() or '-r' in argv:
        call(f'trimmomatic PE -phred33 {input_prefix}_1{EXT} {input_prefix}_2{EXT} {output_prefix}_1_paired{EXT} {output_prefix}_1_unpaired{EXT} {output_prefix}_2_paired{EXT} {output_prefix}_2_unpaired{EXT} ILLUMINACLIP:{adapter_path}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36', shell=True)
        log(f'{acc} trimado.', __file__)
    else:
        print (f"'{str(output_sample_path)}' já existe. Pulando '{acc}'.")
