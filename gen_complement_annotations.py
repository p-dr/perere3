# description: invert heads' and genes' strands on feature's gff.
# in: pardir/'genome_annotation/head_annotations.gff3' 
# in: pardir/'genome_annotation/gene_annotations.gff3' 
# out: pardir/'genome_annotation/head_complement_annotations.gff3' 
# out: pardir/'genome_annotation/gene_complement_annotations.gff3' 

import pandas as pd
from utils import redo_flag, pardir, GFF3_COLUMNS, safe_open

def main():
    for kind, pattern in (('head', r'head\d+'), ('gene', r'Smp_\d+')):
        out_path = pardir/f'genome_annotation/{kind}_complement_annotations.gff3'
        outfile = safe_open(out_path, exist_ok=False)
        gff = pd.read_table(pardir/f'genome_annotation/{kind}_annotations.gff3', names=GFF3_COLUMNS)

        print(gff.strand.head())
        gff.loc[gff.strand == '+', 'strand'] = 'plus'
        gff.loc[gff.strand == '-', 'strand'] = '+'
        gff.loc[gff.strand == 'plus', 'strand'] = '-'
        print(gff.strand.head())

        gff['attributes'] = gff.attributes.str.replace(
            pattern,
            lambda match: match.group(0) + '_complement',
            regex=True)

        gff.to_csv(outfile, sep='\t', header=False, index=False)
        print(f"Wrote to '{str(out_path)}'.")
        outfile.close()
 
if __name__ == '__main__':
    main()
