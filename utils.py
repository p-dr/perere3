from pathlib import Path

scripts_path = Path(__file__).parent
pardir = scripts_path.parent
genome_path = pardir/'seqs/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa'

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

try:
    main_name = Path(__main__.__file__).stem
except AttributeError:
   main_name = '(shell)' 

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
    # WARNING: if attribute has multiple values (eg. 'key=value1; value2; value3;'),
    # only first value is considered.
    expanded_attr = expanded_attr.apply(lambda l: dict([i.split('=') for i in l if '=' in i]))
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


def safe_open(path, mode='w', exist_ok=True):
    message = f"'{path}' já existe, nada será feito. Use '-r' se quiser sobrescrever."

    if path.exists() and not redo_flag:
        if exist_ok:
            print(message)
        else:
            raise FileExistsError(message)

    else:
        return open(path, mode)


import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('GTK3Agg')

grafdir = pardir/'graficos'
oldgrafdir = grafdir/'old'
for folder in (grafdir, oldgrafdir):
    folder.mkdir(exist_ok=True)
figcount = 0


def save_all_figs():
    global figcount
    timestamp = dt.now().strftime('%Y-%-m-%-d-%Hh%Mm%Ss')

    for fignum in plt.get_fignums():
        plt.figure(fignum)
        save_kwargs = dict(bbox_inches='tight', # pad_inches=0,
                           dpi=300)

        plt.savefig(str(oldgrafdir/f'{main_name}_{figcount}_{timestamp}.png'),
                    **save_kwargs)
        plt.savefig(str(grafdir/f'{main_name}_{figcount}.png'),
                    **save_kwargs)
        figcount += 1


def plot_box(count_df, x, y, bins=20, *args):
    newy = str(y)+'_bins'
    count_df.loc[:, newy] = pd.cut(count_df[y], bins)
    count_df.boxplot(x, newy, rot=90, showfliers=False, *args)
    plt.tight_layout()


def boxplot(data):
    """data is list of list of values"""
    for j, d in enumerate(data):
        c = 'C' + str(j)
        plt.boxplot(d, widths=.75, positions=[j], showfliers=False,
                    patch_artist=True,
                    boxprops=dict(facecolor=c, color=c),
                    capprops=dict(color=c),
                    whiskerprops=dict(color=c),
                    flierprops=dict(color=c, markeredgecolor=c),
                    medianprops=dict(solid_capstyle='projecting', color='black')
                    )


from scipy.stats import mannwhitneyu

def box_compare(a, b, labels=None):
    if labels is None:
        labels = a.name, b.name
    pvalue = mannwhitneyu(a, b).pvalue
    plabel = f'p-valor:\n{pvalue:.5}'
    title = ' / '.join(labels)
    median_ratio = a.median() / b.median()
    print(f'{title:<40}', '|  p-value:',
          pvalue, ['x', 'o'][pvalue < .05],
          '| Median ratio:', median_ratio, sep='\t')

    boxplot([a, b])

    plt.annotate(plabel, (.5, .885), xycoords='axes fraction', ha='center')
    labels = [l + '\n' + str(len(data)) for l, data in zip(labels, [a, b])]
    plt.xticks([0, 1], labels=labels)
    
    return pvalue


def get_subsets(d, cols=None):
    d = d.dropna(subset=['neighbor_gene'])  # discard heads with no NG

    if cols is None:
        cols = d.columns

    subsets = dict()
    subsets['-->->'] = d.loc[(d.gene_stream == 'gh') & d.same_strand, cols]
    subsets['<--<-'] = d.loc[(d.gene_stream == 'hg') & d.same_strand, cols]
    subsets['--><-'] = d.loc[(d.gene_stream == 'gh') & ~d.same_strand, cols]
    subsets['<--->'] = d.loc[(d.gene_stream == 'hg') & ~d.same_strand, cols]
    return subsets
