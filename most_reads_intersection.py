# description: (CLARIFICATION NEEDED, UNCERTAIN AFFIRMATION) Plots something to test if heads with most reads paired to it appears consintently across SRA libraries.
# in: pardir/'genome_annotation/heads_motherlength.tsv'
# out: 
# plot: 

from pandas import read_csv, DataFrame, cut, pivot_table
from matplotlib import pyplot as plt
from utils import pardir
from labutils import rgb_to_hex
from pathlib import Path
from random import random as rd, seed
from numpy import log2
from sys import argv


def add_motherlenghts(count_df):
    motherlengths = read_csv(str(pardir/'genome_annotation/heads_motherlength.tsv'), sep='\t', header=None, names=['head_id', 'motherlength'])
    count_df = count_df.join(motherlengths.set_index('head_id'), on='feature')
    return count_df


def color_from_id(name):

    if name.startswith('head'):
        num = int(name.strip('head'))

    elif name.startswith('_'):
        num=0

    else: # name is a gene name.
        num = int(name.split('_')[1])

    seed(num)
    return rgb_to_hex(int(rd()*255) for i in range(3))


def get_sum(df):
    # use empirical indexes instead of boolean filters as below.
    df = df.loc[df['feature'] != '__not_aligned']
    sum_df = int(df['count'].sum())
    df = df.loc[[not i.startswith('_') for i in df['feature']]]
    return sum_df, df


def plot_hist(count_df, count_mode):
    count_df.hist('log2_count', bins=100, density=True, figsize=(16, 9))
    plt.title(stem.replace('_', ' ')+'s')
    plt.xlabel('$\\log_2$(Contagem de reads + 1)')
    plt.ylabel('Frequência relativa')
    plt.show()
    

def plot_bar(count_df, count_mode, stem):
    N_FIRST = 40
    count_df = count_df.sort_values(count_mode+'_count', ascending=False)[:N_FIRST]

    colors = count_df['feature'].apply(color_from_id).values
    count_df.plot('feature', count_mode+'_count', 'bar',
                  color=colors, title=stem,
                  legend=False, figsize=(16, 9))
    plt.tight_layout()
    plt.show()


def plot_line(count_df, count_mode, stem):
    N_FIRST = 100
    count_df = count_df.sort_values(count_mode+'_count', ascending=False)[:N_FIRST]
    count_df.plot('feature', count_mode+'_count', title=stem, legend=False, figsize=(16, 9))
    plt.tight_layout()
    plt.show()

    
def plot_motherlength(count_df, count_mode, stem):
    count_df = add_motherlenghts(count_df)
    colors = count_df['feature'].apply(color_from_id).values

    count_df.plot('motherlength', count_mode+'_count', 'scatter', color=colors, title=stem, figsize=(16, 9), alpha=.5)
    plt.show()


def plot_box(count_df, count_mode, stem):
    BINS = 50
    count_df = add_motherlenghts(count_df)
    #count_df = count_df[count_df['count'] > 10]
    
    count_df.dropna(inplace=True)
    count_df['ml_bins'] = cut(count_df['motherlength'], BINS)
    count_df.boxplot('count', 'ml_bins')
    plt.title(stem)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

    
def main_switch(count_df, stem):

    print('===============', stem)
    count_df = add_mode_cols(count_df)
    BINS = 100
    MIN_COUNT = 2
    
    count_df = add_motherlenghts(count_df)
    count_df['ml_bins'] = cut(count_df['motherlength'], BINS)
    most_reads = count_df[count_df['ml_bins'] == count_df['ml_bins'][0]]
    most_reads = most_reads[most_reads['count'] >= MIN_COUNT]

    for head in most_reads['feature']:
        file_to_head[head] = file_to_head.get(head, []) + [stem]

        

def add_mode_cols(count_df):
    # Todo: Only add necessary modes.
    count_sum, count_df = get_sum(count_df)
    count_df['abs_count'] = count_df['count']
    count_df['log2_count'] = log2(count_df['count']+1)
    count_df['permillion_count'] = count_df['count'] * 1e6 / count_sum
    count_df['log2permillion_count'] = log2(count_df['permillion_count']+1)
    return count_df


#======================== MAIN ========================#

counted_reads_dir = pardir/'counted_reads'
file_to_head = {}
total_count_heads = DataFrame()

BINS = 100
MIN_COUNT = 2
n_files = 0

for count_file in counted_reads_dir.glob('*'):
    stem = Path(count_file).stem
    count_df = read_csv(count_file, sep='\t', header=None,
                        names=['feature', 'count'])

    if stem.endswith('head'):
        total_count_heads = total_count_heads.add(count_df.set_index('feature'), fill_value=0)
        count_df = add_mode_cols(count_df)

        count_df = add_motherlenghts(count_df)
        count_df['ml_bins'] = cut(count_df['motherlength'], BINS)
        most_reads = count_df[count_df['ml_bins'] == count_df['ml_bins'][0]]
        most_reads = most_reads[most_reads['count'] >= MIN_COUNT]

        for head in most_reads['feature']:
            file_to_head[head] = file_to_head.get(head, []) + [stem]

        n_files += 1


file_count = [(i, len(j)/n_files*100) for i, j in file_to_head.items()]
file_count = DataFrame(file_count, columns=['Head ID', 'Porcentagem de arquivos em que aparece'])
file_count.sort_values(file_count.columns[1], ascending=False, inplace=True)

print(file_count)
file_count.plot(*file_count.columns, 'bar', legend=False, title=file_count.columns[1]+', total de arquivos: '+str(n_files))
plt.tight_layout()
plt.show()
        

