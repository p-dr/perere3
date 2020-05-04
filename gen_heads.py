# description: Gets HEAD_LEN posterior bases in S. mansoni genome for each perere3 to genome alignment, generating now called head sequences. Annotates lenght of the parent alignment for each head (motherlength).
# in: inpath u.genome_path
# out: heads_annotations_path heads_outpath motherlength_path
# plot:

from Bio.SeqIO import to_dict, parse
import matplotlib.pyplot as plt
from pandas import read_table, DataFrame
from sys import argv
import utils as u
from sr3_filter import filtered_outpath as inpath

# A filtragem por evalue é feita no próprio BLAST.
# How much does the GTAA repeats influnces the final result?

HEAD_LEN = 200  # final head length
GTAA_WINDOW_LEN = 2000  # search window length for GTAA repetitions
PREFIX_LEN = 12  # Length of printed transposon tail in verbose mode
# Maximum allowed difference between total known length of Perere-3 and the length of the alignment being considered. Intends to discard heads of truncated copies, as truncation tends to occur at 5' end, removing the copy promoter.
MAX_DISTANCE_FROM_END = 5

heads_annotations_path = u.pardir/'genome_annotation/head_annotations.gff3'
heads_outpath = u.pardir/'seqs/heads.fa'
motherlength_path = u.pardir/'genome_annotation/heads_motherlength.tsv'


def find_gtaa_break(seq, phase_tol=3, *args, **kwargs):
    """phase_tol is the phase tolerance of analysis. Warning: it multiplies processing time."""
    return max(phase + gtaa_break_pos(seq[phase:],*args, **kwargs)
               for phase in range(phase_tol + 1))


def gtaa_break_pos(seq, max_dirty_blocks=2, max_errors=2, verbose=False):
    """
    Returns position of seq string where 'GTAA' repetitiveness is broken.
    Tolerates at most max_errors in max_dirty_blocks blocks of four characters.
    """
    seq = seq.upper()
    errors = 0
    count_dirty_blocks = 0

    for i in range(0, len(seq), 4):

        for bp1, bp2 in zip(seq[i:i+4], 'GTAA'):
            if bp1 != bp2:
                errors += 1
                if verbose:
                    print (bp1.lower(), end='')

            elif verbose:
                print(bp1, end='')

            if errors > max_errors:
                ret = i-count_dirty_blocks*4
                if verbose:
                    print ('\n'+seq[:ret]+'|'+seq[ret:])
                return ret

        if errors:
            count_dirty_blocks += 1

            if count_dirty_blocks == max_dirty_blocks:
                count_dirty_blocks = 0
                errors = 0

    if verbose:
        print('\nNo break position found. Returning len(seq).')
    return len(seq)



def main():
    truncated_count = 0
    #======================== LEITURA ========================#
    heads_annotations_file = u.safe_open(heads_annotations_path, exist_ok=False)
    heads_outfile = u.safe_open(heads_outpath, exist_ok=False)
    motherlength_outfile = u.safe_open(motherlength_path, exist_ok=False)

    print('Lendo alinhamentos filtrados do Perere3...', end=' ')
    filtered_perere3_vs_genoma = read_table(inpath)
    print(f"'{inpath}' lido.")

    print('Lendo genoma de S. mansoni...', end=' ')
    genomedict = to_dict(parse(str(u.genome_path), 'fasta'))
    print('Dicionário criado.')



    #======================== GET HEADS ========================#
    print('Buscando as sequências head no genoma...', end=' ')

    with (u.pardir/'seqs/perere3.fa').open() as per_file:
        perere_len = len(''.join([l.strip() for l in per_file.readlines()][1:]))

    heads = []

    for index, row in filtered_perere3_vs_genoma.iterrows():

        # Discard copies without 3' end.
        if abs(row['qend'] - perere_len) < MAX_DISTANCE_FROM_END:
            genome_piece = genomedict[row['saccver']].seq
            plus_sense = row['sstart'] < row['send']

            if plus_sense:
                head_slice = slice(row['send'], row['send'] + GTAA_WINDOW_LEN)
                proto_head = genome_piece[head_slice]

                if u.verbose:
                    prefix = genome_piece[head_slice.start-PREFIX_LEN : head_slice.start]

            else:
                head_slice = slice(row['send'] - GTAA_WINDOW_LEN - 1, row['send'] - 1)
                proto_head = genome_piece[head_slice].reverse_complement()

                if u.verbose:
                    prefix = genome_piece[head_slice.stop : head_slice.stop + PREFIX_LEN].reverse_complement()

            ###### Algumas dão xabu (???)
            if head_slice.start < 0 or head_slice.stop < 0:
                u.prinf(f'Head descartada com posições:', head_slice)
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

            if u.verbose:
                print(['-', '+'][plus_sense], prefix, proto_head[:skip_gtaa]+' | '+head[:30-skip_gtaa]+'...',
                  f" {len(head)}bp\t{row['pident']:.2f}%\t{row['evalue']:.2e}\t{row['bitscore']:5}\t{row['saccver']}\t{head_slice.start}-{head_slice.stop}")

        else:
            truncated_count += 1


    u.log('Written:', heads_outpath, heads_annotations_path, sep='\n\t')
    u.log(f'{filtered_perere3_vs_genoma.shape[0]} alignments considered.')
    u.log(truncated_count, 'heads discarded as truncated.')

    heads_annotations_file.close()
    heads_outfile.close()
    motherlength_outfile.close()

    #======================== PLOTAR HISTOGRAMAS ========================#

    if u.plot_flag:
        heads_df = DataFrame.from_dict(dict(enumerate(zip(*heads))))

        for j in range(8):
            plt.figure(figsize=(16,9))
            for i in range(12):
                plt.subplot(3,4,i+1)
                heads_df[i+j*12].value_counts().plot(kind='bar')

            plt.show()


if __name__ == '__main__':
    main()
