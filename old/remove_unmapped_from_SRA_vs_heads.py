# Remove lines with '*' (not aligned) from SAM files.
from utils import pardir, redo_flag, log

indir = pardir/'alinhamentos/SRA_vs_heads'
outdir  = pardir/'alinhamentos/SRA_vs_heads/only_mapped'

for sam_path in indir.glob('*.sam'):

    outpath = outdir/sam_path.name
    
    if outpath.exists() and not redo_flag:
        print(f"'{outpath.name}' já existe, nada será feito. Use '-r' para sobrescrever.")
        continue

    with sam_path.open('r') as sam_file:
        with outpath.open('w') as outfile:
            
            for sam_line in sam_file:
                if '*' not in sam_line and not sam_line.startswith('@'):
                    row = sam_line.split('\t')
                    outfile.write(f'{row[2]}\t{row[3]}\t{row[5]}\n')

            print(f"{sam_path.stem} teve reads que não alinharam removidos.")
            log(f"{sam_path.stem} teve reads que não alinharam removidos.", __file__)
