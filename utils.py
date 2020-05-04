from pathlib import Path

scripts_dir = Path(__file__).resolve().parent
plot_scripts_dir = scripts_dir/'plot_scripts'
pardir = scripts_dir.parent
genome_path = pardir/'seqs/sm_genome.fa'

LOG_DIR = pardir/'logs'
GRAF_DIR = pardir/'graficos'
OLDGRAF_DIR = GRAF_DIR/'old'

for folder in (LOG_DIR, GRAF_DIR, OLDGRAF_DIR):
    folder.mkdir(exist_ok=True)

import __main__

try:
    main_name = Path(__main__.__file__).stem
except AttributeError:
   main_name = '(shell)'

LOG_PATH = (LOG_DIR/main_name).with_suffix('.log')


GFF3_COLUMNS = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
BL_COLUMNS = ['qaccver', 'saccver', 'qstart', 'qend', 'sstart', 'send', 'length', 'pident', 'evalue', 'bitscore']

import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('-r', '--redo', action='store_true', help='overwrite existing files')
argparser.add_argument('-v', '--verbose', action='store_true', help='increase verbosity')
argparser.add_argument('-s', '--show', action='store_true', help='show instead of save created plots')
argparser.add_argument('-ss', '--show_save', action='store_true', help='show and save created plots')
argparser.add_argument('-p', '--plot', action='store_true', help='show plot information')
argparser.add_argument('-c', '--clean', action='store_true', help='delete large files after processing them')
args = argparser.parse_args()

# only for back compatibility
redo_flag = args.redo
verbose = args.verbose
show_flag = args.show
plot_flag = args.plot


def prinf(text, *args, **kwargs):
    if verbose:
        print(text, *args, **kwargs)


from datetime import datetime as dt

def log(*args, tee=True, **kwargs):
    """Writes text to log file in pardir/logs.
    """
    date = dt.now().strftime('%c')

    with LOG_PATH.open('a') as log_file:
        print(f"{date}:", *args, file=log_file, **kwargs)
        if tee:
            print(*args, **kwargs)


def overlaps(*args):
    """ Test if two closed intervals (boundaries included) respresented
    by two tuples of the form '(start pos., end pos.)' overlap. As this
    is designed to deal with DNA sequences, it can only return True if
    boundary values are on the same sense.
    """
    if len(args) == 1:
        args = args[0]
    if len(args) == 2:
        t1, t2 = args
    elif len(args) == 4:
        t1, t2 = args[:2], args[2:]

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


from time import time

# timeit does exactly the same, but better.
def time_func(func, nt=1000, *args, **kwargs):
    t0 = time()
    for i in range(nt):
        func(*args, **kwargs)
    return time()-t0


from pandas import DataFrame
import pandas as pd

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
            if exist_ok == 'exit':
                exit()
        else:
            raise FileExistsError(message)

    else:
        return open(path, mode)


import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm
#matplotlib.use('GTK3Agg')

figcount = 0


def save_all_figs():
    global figcount

    if args.show or args.show_save:
        plt.show()
        if args.show:
            return
        
    timestamp = dt.now().strftime('%Y-%-m-%-d-%Hh%Mm%Ss')

    for fignum in tqdm(plt.get_fignums(), desc='Saving figures'):
        plt.figure(fignum)
        save_kwargs = dict(bbox_inches='tight', # pad_inches=0,
                           dpi=300)

        plt.savefig(str(OLDGRAF_DIR/f'{main_name}_{figcount}_{timestamp}.png'),
                    **save_kwargs)
        plt.savefig(str(GRAF_DIR/f'{main_name}_{figcount}.png'),
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
    label_string = [l + '\n' + str(len(data)) for l, data in zip(labels, [a, b])]
    plt.xticks([0, 1], labels=label_string)
    
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


def multibox_compare(ds, labels=None, arch_height=None, margin=None):
    if labels is None:
        labels = [d.name for d in ds]
    labels = [f'{l}\n{len(d)}'for l, d in zip(labels, ds)]
    print('labels:', labels)
    print('lengths from multibox:', *[len(d) for d in ds])

    boxplot(ds)
    plt.xticks(range(len(ds)), labels=labels)
    height = plt.axis()[-1]
    fig_height = height - plt.axis()[-2]

    if arch_height is None:
        arch_height = .02 * fig_height

    if margin is None:
        margin = .04 * fig_height

    for i in range(len(ds)):
        for j in range(i+1, len(ds)):
            plt.plot((i, i, j, j),
                    (height,
                     height + arch_height,
                     height + arch_height,
                     height), '-k', lw=1)

            pvalue = mannwhitneyu(ds[i], ds[j]).pvalue
            plabel = f'p: {pvalue:.5}'
            print(plabel)
            plt.annotate(plabel, ((i+j)/2, height + arch_height),
                         ha='center', va='bottom')
            height += arch_height + margin


import os

def print_header(*args, sep=' ', log=False, **kwargs):
    outs = [LOG_PATH.open('a'), None][not log:]  # You let an open file.
    cols, _ = os.get_terminal_size()
    args = [a.upper() if type(a) == str else str(a) for a in args]
    argstr = f'{sep.join(args):^{cols}}'

    for outfile in outs:
        print('=' * cols, file=outfile)
        print(argstr, file=outfile, **kwargs)
        print('=' * cols, file=outfile)
        print(file=outfile)


def clean(*paths):
    if args.clean:
       [path.unlink() for path in paths] 
       log('Deleted following files:', *paths, sep='\n\t')


def ask_to_call(func, desc, *args, **kwargs):
    # implement not_func to call if inp in 'nN'
    while True:
        inp = input('Delete leftover files? [y/n]:')
        if inp in 'yY':
            func(*args, **kwargs)
            return True
        elif inp in 'nN':
            return False


def user_confirms(desc='Confirm action?', yes_options='yY', no_options='nN'):
    while True:
        inp = input(f'{desc} [{yes_options}/no_options]:')
        if inp in yes_options:
            return True
        elif inp in no_options:
            return False
