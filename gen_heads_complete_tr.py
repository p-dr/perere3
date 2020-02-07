# description: Gets HEAD_LEN posterior bases in S. mansoni genome for each perere3 to genome alignment, generating now called head sequences. Annotates lenght of the parent alignment for each head (motherlength).
# in: pardir/'alinhamentos/filtered_perere3complete_vs_genoma.bl' pardir/'seqs/sm_genome.fa'
# out: pardir/'seqs/heads.fa' pardir/'genome_annotation/head_annotations.gff3' pardir/'genome_annotation/heads_motherlength.tsv'
# plot:

from Bio.SeqIO import to_dict, parse
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
from sys import argv
from utils import find_gtaa_break, pardir, verbose, prinf, genome_path

# A filtragem por evalue é feita no próprio BLAST.

HEAD_LEN = 100                # final head length
GTAA_WINDOW_LEN = 500         # search window length for GTAA repetitions
PREFIX_LEN = 12

if __name__ == '__main__':
    #======================== LEITURA ========================#

    print('Lendo alinhamentos filtrados do Perere3...')
    inpath = pardir/'alinhamentos/filtered_perere3complete_vs_genoma.bl'
    filtered_perere3_vs_genoma = read_csv(inpath, sep='\t')
    print(f"'{inpath}' lido.")

    print('Lendo genoma de S. mansoni...')
    genomedict = to_dict(parse(str(genome_path), 'fasta'))
    print('Dicionário criado.')


    print ('Abrindo arquivos de output e anotação...')

    heads_annotations_path = pardir/'genome_annotation/head_annotations.gff3'
    heads_outpath = pardir/'seqs/heads.fa'
    motherlength_path = pardir/'genome_annotation/heads_motherlength.tsv'

    heads_annotations_file = heads_annotations_path.open('w')
    # heads_annotations_file.write('\t'.join(GFF3_COLUMNS)+'\n') # Write header.
    
    heads_outfile = heads_outpath.open('w')
    motherlength_outfile = motherlength_path.open('w')
    
    #======================== GET HEADS ========================#

    print('Buscando as sequências head no genoma...')
    
    with (pardir/'seqs'/'perere3complete.fa').open() as per_file:
        perere_len = len(''.join([l.strip() for l in per_file.readlines()][1:]))

    MAX_DISTANCE_FROM_END = 3
    heads = []

    for index, row in filtered_perere3_vs_genoma.iterrows():

        if abs(row['qend'] - perere_len) < MAX_DISTANCE_FROM_END:
            genome_piece = genomedict[row['saccver']].seq
            plus_sense = row['sstart'] < row['send']

            if plus_sense:
                head_slice = slice(row['send'], row['send'] + GTAA_WINDOW_LEN)
                proto_head = genome_piece[head_slice]
                
                # verbose use only
                prefix = genome_piece[head_slice.start-PREFIX_LEN : head_slice.start]

            else:
                head_slice = slice(row['send'] - GTAA_WINDOW_LEN - 1, row['send'] - 1)
                proto_head = genome_piece[head_slice].reverse_complement()

                # verbose use only
                prefix = genome_piece[head_slice.stop : head_slice.stop + PREFIX_LEN].reverse_complement()

            ###### Algumas dão xabu
            if head_slice.start < 0 or head_slice.stop < 0:
                prinf(f'Head descartada com posições:', head_slice)
                continue
            
            skip_gtaa = find_gtaa_break(proto_head)
            head = proto_head[skip_gtaa : skip_gtaa + HEAD_LEN]


            #======================== ANOTAR HEAD NO GFF ========================#

            if plus_sense:
                start_pos = row['send'] + 1 + skip_gtaa
                end_pos = start_pos + HEAD_LEN

            else:
                end_pos = row['send'] - 1 - skip_gtaa
                start_pos = end_pos - HEAD_LEN

            heads_annotations_file.write('\t'.join([row['saccver'], 'WormBase_imported', 'gene',
                                                    str(start_pos), str(end_pos), '.', ['-', '+'][plus_sense], '.',
                                                    f'gene_id=head{index};motherlength={row["length"]};length={HEAD_LEN}'])+'\n')

            motherlength_outfile.write(f"head{index}\t{row['length']}\n")

            #======================== ESCREVER HEAD EM FASTA ========================#

            heads_outfile.write(f'>head{index}\n' + str(head) + '\n')

            #========================================================================#


            prinf(['-', '+'][plus_sense], prefix, proto_head[:skip_gtaa]+' | '+head[:30-skip_gtaa]+'...',
                  f" {len(head)}bp\t{row['pident']:.2f}%\t{row['evalue']:.2e}\t{row['bitscore']:5}\t{row['saccver']}\t{head_slice.start}-{head_slice.stop}")


    print('Todas as sequências encontradas.')
    print(f"Heads salvas em '{str(heads_outpath)}' e anotadas em '{str(heads_annotations_path)}'.")
    print ('Fechando arquivos e encerrando.')

    heads_annotations_file.close()
    heads_outfile.close()

    #======================== PLOTAR HISTOGRAMAS ========================#

    if '--plot' in argv:

        heads_df = DataFrame.from_dict(dict(enumerate(zip(*heads))))

        for j in range(8):
            plt.figure(figsize=(16,9))
            for i in range(12):
                plt.subplot(3,4,i+1)
                heads_df[i+j*12].value_counts().plot(kind='bar')

            plt.show()
