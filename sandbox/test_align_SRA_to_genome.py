import sys
sys.path.append("..")

import align_SRA_to_genome as align
import count_SRA_reads as count
import subprocess as sp
from pathlib import Path
import utils as u

print('WARNING: Always run from pardir/test directory.')
logs_dir = Path('./logs')
alignments_dir = Path('./alignments')
trimmed_dir = Path('./trimmed')
counts_dir = Path('./counts')

for d in logs_dir, alignments_dir, trimmed_dir, counts_dir:
    d.mkdir(exist_ok=True)

u.LOG_DIR = logs_dir
align.out_dir = alignments_dir
align.trimmed_data_dir = trimmed_dir
count.out_dir = counts_dir

acc = 'SRX7450678_head100000'

if len(sys.argv) < 2:
    print("'a' para alinhar e 'c' contar.")
    exit()

#if 'a' in sys.argv[-1]:
if u.redo_flag:
    align.align_acc(acc)

#if 'c' in sys.argv[-1]:
if u.verbose:
    sam = list(align.out_dir.glob(acc + '*'))[0]
    print(f'Contando reads do arquivo {str(sam)}.')
    count.count_reads((sam, 'gene'))
    count.count_reads((sam, 'gene_complement'))

