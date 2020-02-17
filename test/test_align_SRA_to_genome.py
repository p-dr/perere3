import sys
sys.path.append("..")

import align_SRA_to_genome as align
import count_SRA_reads as count
import subprocess as sp
from pathlib import Path
import utils as u

align.out_dir = Path('./alignments')
count.out_dir = Path('./counts')

acc = 'ERR022874'

if len(sys.argv) < 2:
    print("'a' para alinhar e 'c' contar.")
    exit()

if 'a' in sys.argv[-1]:
    #align.align_acc(acc)
    sp.run(['hisat2',
        '-x', u.pardir/'sm_genome_ht2',
        '--rna-strandness', 'FR',
        '--sra-acc', acc,
        '-S', f'./alignments/{acc}_online.sam',
    ])

sam = list(align.out_dir.glob('*'))[0]

if 'c' in sys.argv[-1]:
    count.count_reads((sam, 'gene'))
    count.count_reads((sam, 'gene_complement'))

