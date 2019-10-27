from pathlib import Path

scripts_path = Path(__file__).parent
pardir = scripts_path.parent

GFF3_COLUMNS = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
BL_COLUMNS = ['qaccver', 'saccver', 'qstart', 'qend', 'sstart', 'send', 'length', 'pident', 'evalue', 'bitscore']

from sys import argv

redo_flag = '-r' in argv
if redo_flag:
    argv.remove('-r')

verbose = '-v' in argv
if verbose:
    argv.remove('-v')

show_flag = '-s' in argv
if show_flag:
    argv.remove('-s')

progress_flag = '-p' in argv
if progress_flag:
    argv.remove('-p')

def prinf(text, *args, **kwargs):
    global verbose
    if verbose:
        print(text, *args, **kwargs)

from datetime import datetime as dt
import __main__
main_name = Path(__main__.__file__).stem

def log(text, author_script=main_name):
    """Writes text to log file in pardir/logs.
    """
    log_path = pardir/'logs'/(author_script+'.log')
    date = dt.now().strftime('%c')

    with log_path.open('a') as log_file:
        log_file.write(f"{date}: {text}\n")


def overlaps(t1, t2):
    """ Test if two closed intervals (boundaries included) respresented
    by two tuples of the form '(start pos., end pos.)' overlap. As this
    is designed to deal with DNA sequences, it can only return True if
    boundary values are on the same sense.
    """
    if isinstance(t1, int):
        t1 = (t1, t1)

    if isinstance(t2, int):
        t2 = (t2, t2)

    t1ascending = t1[1] >= t1[0]
    t2ascending = t2[1] >= t2[0]

    # test if both are on the same strand (both start/end
    # positions are either ascending or descending)
    if t1ascending == t2ascending:

        t1 = sorted(t1)
        t2 = sorted(t2)

        return not ((t1[0] > t2[1]) or (t1[1] < t2[0]))

    # if not, comparison makes no sense
    else:
        return False


def find_gtaa_break(seq, max_dirty_blocks=2, max_errors=2, verbose=False):
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


from time import time

def time_func(func, nt=1000, *args, **kwargs):
    t0 = time()
    for i in range(nt):
        func(*args, **kwargs)
    return time()-t0

from pandas import read_csv, DataFrame
import pandas as pd

def read_tsv(*args, **kwargs):
    return read_csv(*args, sep='\t', **kwargs)


def parse_gff_attributes(attr, id_as_index=True, gene_id='gene_id'):
    expanded_attr = attr.str.split(';')
    expanded_attr = expanded_attr.apply(lambda l: dict([i.split('=') for i in l]))
    ret = DataFrame(list(expanded_attr))
    if id_as_index:
        ret.set_index(gene_id, inplace=True)
    return ret.infer_objects()


def unfold_gff(df):
    attr = df.attributes
    df = df.drop('attributes', 1)
    attr = parse_gff_attributes(attr, id_as_index=False)
    ret = pd.concat([df, attr], 1)
    return ret


def safe_open(path, mode):
    if path.exists() and not redo_flag:
        print(f"'{path}' já existe, nada será feito. Use '-r' se quiser sobrescrever.")
        exit()

    return open(path, mode)


import matplotlib.pyplot as plt
grafdir = pardir/'graficos'
figcount = 0


def save_all_figs():
    global figcount
    timestamp = dt.now().strftime('%Y-%-m-%-d-%Hh%Mm%Ss')

    for fignum in plt.get_fignums():
        plt.figure(fignum, dpi=300)
        plt.savefig(grafdir/f'old/{main_name}_{figcount}_{timestamp}.png')
        plt.savefig(grafdir/f'{main_name}_{figcount}.png')
        figcount += 1


def plot_box(count_df, x, y, bins=20, *args):
    newy = str(y)+'_bins'
    count_df[newy] = pd.cut(count_df[y], bins)
    count_df.boxplot(x, newy, rot=90, showfliers=False, *args)
    plt.tight_layout()
