# description: Plots how many reads aligned to each head, under many plotting formats (--help).
# in: pardir/'genome_annotation/heads_motherlength.tsv' pardir/'counted_reads'
# out: 
# plot: 

from pandas import read_csv, DataFrame, cut, pivot_table
from matplotlib import pyplot as plt
from utils import pardir, save_all_figs
from labutils import rgb_to_hex
from pathlib import Path
from random import random as rd, seed
import numpy as np
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
    

def plot_bar(count_df, count_mode, stem):
    N_FIRST = 40
    count_df = count_df.sort_values(count_mode+'_count', ascending=False)[:N_FIRST]

    colors = count_df['feature'].apply(color_from_id).values
    count_df.plot('feature', count_mode+'_count', 'bar',
                  color=colors, title=stem,
                  legend=False, figsize=(16, 9))
    plt.tight_layout()


def plot_line(count_df, count_mode, stem):
    N_FIRST = 100
    count_df = count_df.sort_values(count_mode+'_count', ascending=False)[:N_FIRST]
    count_df.plot('feature', count_mode+'_count', title=stem, legend=False, figsize=(16, 9))
    plt.tight_layout()
    
    
def plot_motherlength(count_df, count_mode, stem):
    count_df = add_motherlenghts(count_df)
    colors = count_df['feature'].apply(color_from_id).values

    count_df.plot('motherlength', count_mode+'_count', 'scatter', color=colors,
                  title=stem, figsize=(16, 9), alpha=.5, figure=30)


def plot_box(count_df, count_mode, stem):
    BINS = 20
    count_df = add_motherlenghts(count_df)
    # count_df = count_df[count_df['count'] > 10]
    
    count_df['ml_bins'] = cut(count_df['motherlength'], BINS)
    count_df.boxplot('count', 'ml_bins', rot=90, showfliers=False)
    plt.title(stem)
    plt.tight_layout()


def plot_hexbin(count_df, count_mode, stem):
    count_df = add_motherlenghts(count_df)
#    count_df = count_df[count_df['count'] > 10]
    count_df.plot.hexbin('motherlength', count_mode+'_count',
                         gridsize=10, cmap='viridis')


def plot_heatmap(count_df, count_mode, stem, bins = 10):
    count_df = add_motherlenghts(count_df)
    count_df['ml_bins'] = cut(count_df['motherlength'], bins)
    count_df['count_bins'] = cut(count_df[count_mode+'_count'], bins)
    ptab = pivot_table(count_df, values=count_mode+'_count',
                       index='count_bins', columns='ml_bins',
                       aggfunc=np.sum)
    ptab.fillna(0, inplace=True)
    plt.pcolor(ptab)
    plt.colorbar()
    plt.xlabel('motherlength')
    plt.ylabel(count_mode)
    plt.title(stem)

def main_switch(count_df, stem, count_mode):

    count_df = add_mode_cols(count_df)
    
    if 'hist' in argv:
        plot_hist(count_df, count_mode, stem)

    if 'bar' in argv:
        plot_bar(count_df, count_mode, stem)
    
    if 'line' in argv:
        plot_line(count_df, count_mode, stem)

    if stem.endswith('head'):
        
        if ('ml' in argv or 'motherlength' in argv):
            plot_motherlength(count_df, count_mode, stem)

        if 'box' in argv:
            plot_box(count_df, count_mode, stem)

        if 'hex' in argv:
            plot_hexbin(count_df, count_mode, stem)

        if 'heatmap' in argv:
            plot_heatmap(count_df, count_mode, stem)

    plt.show()


def add_mode_cols(count_df):
#Only add necessary modes.
    count_sum, count_df = get_sum(count_df)
    count_df['abs_count'] = count_df['count']
    count_df['log2_count'] = np.log2(count_df['count']+1)
    count_df['permillion_count'] = count_df['count'] * 1e6 / count_sum
    count_df['log2permillion_count'] = np.log2(count_df['permillion_count']+1)
    return count_df


#======================== MAIN ========================#
if __name__ == '__main__':
    
    counted_reads_dir = pardir/'counted_reads'
    count_modes = [mode.strip('-') for mode in argv if mode.startswith('-')]

    if '--help' in argv:
        print("\nModos de formatar a contagem:")
        print("\t-abs\t\tContagem absoluta de reads;")
        print("\t-log2\t\tLog na base 2 da contagem de reads;")
        print("\t-permillion\tContagem por milhão de reads alinhados no total;")
        print("\t-log2permillion\tLog2 da contagem por milhão de reads alinhados no total.")
        print("\nOpções de gráfico:\n\thist;\n\tbar;\n\tline;\n\tml ou motherlength;\n\tbox\t(boxplot);\n\thex\t(hexbin).")
        print("\nSoma:\n\tsum\t\texibe um gráfico ao final com todas as contagens somadas;\n\tonly\t\texibe somente a contagem final total.")
        print("\nAjuda:\n\t--help\t\tsomente imprime esse texto.")
        exit()

    for count_mode in count_modes:

        total_count_heads = DataFrame()
        total_count_genes = DataFrame()

        for count_file in counted_reads_dir.glob('*'):
            stem = Path(count_file).stem

            # ALERTA DE GAMBIARRA
            if stem.startswith('aggregated'):
                continue

            count_df = read_csv(count_file, sep='\t', header=None,
                                names=['feature', 'count'])

            if stem.endswith('gene'):
                total_count_genes = total_count_genes.add(count_df.set_index('feature'), fill_value=0)
            else:
                total_count_heads = total_count_heads.add(count_df.set_index('feature'), fill_value=0)

            if 'only' not in argv:
                main_switch(count_df, stem, count_mode)

        if 'sum' in argv:
            stem_list = ['Contagem total: gene', 'Contagem total: head']
            for i, total_count in enumerate([total_count_genes, total_count_heads]):
                main_switch(total_count.reset_index(), stem_list[i], count_mode)



