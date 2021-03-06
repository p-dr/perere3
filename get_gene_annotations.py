# description: Removes non-gene annotations from S. mansoni genome annotations. Calculates gene length.
# in: pardir/'genome_annotation/sm_annotations.gff3'
# out: pardir/'genome_annotation/gene_annotations.gff3'

from pandas import read_csv
from utils import pardir, redo_flag, GFF3_COLUMNS, unfold_gff
import utils as u
from wget import download
from subprocess import call
import os

raw_annotations_path = pardir/'genome_annotation/sm_annotations.gff3'
outpath = pardir/'genome_annotation/gene_annotations.gff3'

#============== REMOVER ANOTAÇÕES NÃO-GÊNICAS ===============#
def main():
    with u.safe_open(outpath, exist_ok=False) as outfile:
        raw_annotations = read_csv(raw_annotations_path,
                                   sep='\t', comment='#',
                                   header=None, names=GFF3_COLUMNS)


        print('Leitura encerrada. Removendo anotações não-gênicas...')

        genes_gff = raw_annotations.loc[raw_annotations['type'] == 'gene']
        genes_gff.loc[:, ['start', 'end']] = genes_gff[['start', 'end']].astype(int)

        lengths = genes_gff.end - genes_gff.start
        genes_gff.loc[:, 'attributes'] = genes_gff.attributes.str.replace('ID=gene:', 'gene_id=')
        genes_gff['attributes'] = genes_gff.attributes.str.extract(r'(gene_id.*Name[^;]+)')
        # com loc não funciona (?!):
        # genes_gff.loc[:, 'attributes'] = genes_gff.attributes.str.extract(r'(gene_id.*Name[^;]+)')
        genes_gff.loc[:, 'attributes'] += ';length=' + lengths.astype(str)


        # ###### REMOVER GENES COM FIM OU INÍCIO COINCIDENTES
        genes_gff = genes_gff.loc[lengths.sort_values().index]  # Ordenar por tamanho
        # Manter o maior gene entre os que coincidem.
        genes_gff = genes_gff.drop_duplicates(['seqid', 'start'], keep='last')
        genes_gff = genes_gff.drop_duplicates(['seqid', 'end'], keep='last')

        if genes_gff.duplicated(['seqid', 'start']).sum() or genes_gff.duplicated(['seqid', 'end']).sum():
            print('ERRO: HÁ GENES COM INÍCIO/TÉRMINO DUPLICADOS:')
            print(genes_gff[genes_gff.duplicated(['seqid', 'start'], keep=False)])
            print(genes_gff[genes_gff.duplicated(['seqid', 'end'], keep=False)])
            raise ValueError

        genes_gff = genes_gff.sort_values(['seqid', 'start'])
        genes_gff.to_csv(outfile, sep='\t', index=False, header=None)

        print(f"Anotações gênicas mantidas em '{str(outpath)}'.")


if __name__ == '__main__':
    main()
