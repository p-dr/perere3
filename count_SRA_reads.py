# description: Counts how many reads were aligned to each gene and head, using HTSeq. Saves it as an attribute on GFF file.
# in: pardir/'alinhamentos/SRA_vs_genoma' pardir/'genome_annotation/head_annotations.gff3' pardir/'genome_annotation/gene_annotations.gff3'
# out: pardir/'counted_reads'

from glob import iglob
from pathlib import Path
from subprocess import call
from utils import pardir, redo_flag, log

alignment_dir = pardir/'alinhamentos/SRA_vs_genoma'
annotations_dir = pardir/'genome_annotation'
out_dir = pardir/'counted_reads'

for alignment_file in iglob(str(alignment_dir/'*')):
    for kind in ['head', 'gene']:

        acc = Path(alignment_file).stem
        annotations_file = str(annotations_dir/(kind+'_annotations.gff3'))
        out_path = out_dir/(acc+'_'+kind+'.csv')
        
        if not out_path.exists() or redo_flag:
            print (f"Contando reads de {acc} em cada anotação de {kind}...")
            call(f'python -m HTSeq.scripts.count {alignment_file} {annotations_file} -t gene -q > {str(out_path)}', shell=True)
            print (f'\n{out_path.name}: htseq-count terminou de rodar.')
            log(f"'{str(out_path)}' finalizado.", __file__)
        else:
            print (f"'{str(out_path)}' já existe. Pulando '{acc}_{kind}'.")

print ('Pronto.')
