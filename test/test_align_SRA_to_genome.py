import sys
sys.path.append("..")

import align_SRA_to_genome as align
import count_SRA_reads as count
from pathlib import Path

align.out_dir = Path('./alignments')
count.out_dir = Path('./counts')

acc = 'ERR022874'

if len(sys.argv) < 2:
    print("'a' para alinhar e 'c' contar.")
    exit()

if 'a' in sys.argv[-1]:
    align.align_acc(acc)

sam = list(align.out_dir.glob('*'))[0]

if 'c' in sys.argv[-1]:
    count.count_reads((sam, 'gene'))

