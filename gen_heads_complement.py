# description: invert heads' strands on heads' gff.
# in: pardir/'genome_annotation/head_annotations.gff3' 
# out: pardir/'genome_annotation/head_complement_annotations.gff3' 

import pandas as pd
from utils import redo_flag, pardir, GFF3_COLUMNS

out_path = pardir/'genome_annotation/head_complement_annotations.gff3' 
gff = pd.read_table(pardir/'genome_annotation/head_annotations.gff3', names=GFF3_COLUMNS)

print(gff.strand)
gff.loc[gff.strand == '+', 'strand'] = 'plus'
gff.loc[gff.strand == '-', 'strand'] = '+'
gff.loc[gff.strand == 'plus', 'strand'] = '-'
print(gff.strand)

gff.to_csv(out_path, sep='\t', header=False, index=False)
print(f"Wrote to '{str(out_path)}'.")
