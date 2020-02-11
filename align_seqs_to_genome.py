# description: Generates perere3 and SR3 vs genome alignments.
# in: pardir/'seqs/perere3.fa' pardir/'seqs/sr3.fa' pardir/'sm_genome_db'
# out: pardir/'alinhamentos/perere3_vs_genome.bl' pardir/'alinhamentos/sr3_vs_genome.bl'

from subprocess import run
from utils import pardir, safe_open, BL_COLUMNS
from multiprocessing import cpu_count
n_cpu = cpu_count()
# Usar biopython para blastear?

QUERIES = ['perere3', 'sr3']
genomedb_path = pardir/'sm_genome_db/sm_genome_db'

#================== CRIAR ARQUIVOS DE ALINHAMENTO ==================#

for query in QUERIES:

    query_path = pardir/f'seqs/{query}.fa'
    out_path = pardir/f'alinhamentos/{query}_vs_genome.bl'
    out_file = safe_open(out_path)

    if out_file is None:
        continue

    print(f'Procurando alinhamentos de {query} contra genoma...')
    run((f"blastn -task blastn -query {str(query_path)} -db {str(genomedb_path)} "
         f"-outfmt '6 {' '.join(BL_COLUMNS)}' -out {str(out_path)} "
         f"-evalue 1e-10 -num_threads '{n_cpu}'"), shell=True)
    print(f'Alinhamentos salvos em {str(out_path)}.\n')
