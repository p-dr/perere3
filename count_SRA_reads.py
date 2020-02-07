# description: Counts how many reads were aligned to each gene and head, using HTSeq. Saves it as an attribute on GFF file.
# in: pardir/'alinhamentos/SRA_vs_genoma'
# in: pardir/'genome_annotation/head_annotations.gff3'
# in: pardir/'genome_annotation/gene_annotations.gff3'
# in: pardir/'genome_annotation/head_complement_annotations.gff3'
# in: pardir/'genome_annotation/gene_complement_annotations.gff3'
# out: pardir/'counted_reads'

from glob import glob
from pathlib import Path
from subprocess import call
from utils import pardir, redo_flag, log
import multiprocessing as mp
from datetime import datetime

alignment_dir = pardir/'alinhamentos/SRA_vs_genoma_bam'
annotations_dir = pardir/'genome_annotation'
out_dir = pardir/'counted_reads'
out_dir.mkdir(exist_ok=True)

in_format = 'bam'

def count_reads(args):
    alignment_file, kind = args

    acc = Path(alignment_file).stem
    annotations_file = str(annotations_dir/(kind+'_annotations.gff3'))
    out_path = out_dir/(acc+'_'+kind+'.csv')
    
    if not out_path.exists() or redo_flag:
        print (f"Contando reads de {acc} em cada anotação de {kind}...")
        t0 = datetime.now()
        call((f'python -m HTSeq.scripts.count {alignment_file} -f {in_format} '
              f'{annotations_file} -t gene -q > {str(out_path)}'), shell=True)
        dt = datetime.now() - t0
        print (f'\n{out_path.name}: htseq-count terminou de rodar. Tempo decorrido: {dt}')
        log(f"'{str(out_path)}' finalizado em {dt}.", __file__)
    else:
        print (f"'{str(out_path)}' já existe. Pulando '{acc}_{kind}'.")


alignment_files = [f for f in alignment_dir.glob('*.' + in_format) if 'tmp' not in f.name]
print('Arquivos de alinhamento encontrados: ', *alignment_files)
types = ['head', 'gene', 'head_complement', 'gene_complement']

with mp.Pool(processes=mp.cpu_count()) as pool:
    t0 = datetime.now()
    pool.map(count_reads, ((af, t) for af in alignment_files for t in types))

dt = datetime.now() - t0
print('\nTodos os arquivos foram processados.')
print(f'Tempo total decorrido: {dt}')
log('Sessão de contagem encerrada com duração {dt}.')
