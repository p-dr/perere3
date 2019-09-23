# description: Get coordinates of genome's masked low-complexity or repetitive regions.
# in: pardir/'seqs/schistosoma_mansoni.PRJEA36577.WBPS14.genomic_masked.fa'
# out: pardir/'genome_annotation/masked_map.tsv'

from utils import pardir, safe_open, prinf
from Bio.SeqIO import parse
from re import finditer

HEADER = ['seqid', 'start', 'end']

inpath = pardir/'seqs/schistosoma_mansoni.PRJEA36577.WBPS14.genomic_masked.fa'
outpath = pardir/'genome_annotation/mask_map.tsv'

if __name__ == '__main__':

    outfile = safe_open(outpath, 'w')
    outfile.write('\t'.join(HEADER) + '\n')
    genome = parse(inpath, 'fasta')

    for record in genome:
        for match in finditer('N+', str(record.seq)):
            outfile.write('\t'.join([record.id, str(match.start()), str(match.end())]) + '\n')

    outfile.close()
