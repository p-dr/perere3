# description: Generates perere3 and SR3 vs genome alignments.
# in: pardir/'seqs/perere3complete.fa' pardir/'seqs/sr3complete.fa' pardir/'genome_db'
# out: pardir/'alinhamentos/perere3complete_vs_genoma.bl' pardir/'alinhamentos/sr3complete_vs_genoma.bl'

from subprocess import run
from os.path import exists
from utils import pardir, redo_flag, BL_COLUMNS

# Usar biopython para blastear?

# QUERIES = ['perere3', 'sr3', 'perere3complete', 'sr3complete']
QUERIES = ['perere3complete', 'sr3complete']

genomedb_path = pardir/'genome_db/smgenome'

#================== CRIAR ARQUIVOS DE ALINHAMENTO ==================#

for query in QUERIES:
    
    query_path = pardir/f'seqs/{query}.fa'
    out_path = pardir/f'alinhamentos/{query}_vs_genoma.bl'
    out_not_exists = not exists(out_path)
    
    if out_not_exists:
        print(f"'{str(out_path)}' não existe.")
        
    if out_not_exists or redo_flag:
        print(f'Procurando alinhamentos de {query} contra genoma...')
        run(f"blastn -task blastn -query {str(query_path)} -db {str(genomedb_path)} -outfmt '6 {' '.join(BL_COLUMNS)}' -out {str(out_path)} -evalue 1e-10", shell=True)
        print(f'Alinhamentos salvos em {str(out_path)}.\n')
