import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from utils import show_flag

plt.style.use('seaborn')
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#22848e',
                                                    '#f88904',
                                                    '#78c18c'])
mpl.rcParams['savefig.transparent'] = True
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['font.family'] = 'Montserrat'

weight = 'semibold'
mpl.rcParams['font.weight'] = weight
mpl.rcParams['axes.titleweight'] = weight
mpl.rcParams['axes.labelweight'] = weight
mpl.rcParams['figure.titleweight'] = weight

gray = '#e2e2e2'
mpl.rcParams['axes.labelcolor'] = gray
mpl.rcParams['axes.edgecolor'] = gray
mpl.rcParams['text.color'] = gray

brighter = '#6a624a'
mpl.rcParams['xtick.color'] = brighter
mpl.rcParams['ytick.color'] = brighter
mpl.rcParams['grid.color'] = brighter

bgcolor = '#343025'
# mpl.rcParams['boxplot.patchartist'] = True
mpl.rcParams['boxplot.showfliers'] = False
# mpl.rcParams['boxplot.boxprops.color'] = 'C0'
mpl.rcParams['boxplot.boxprops.linewidth'] = 3
mpl.rcParams['boxplot.whiskerprops.linewidth'] = 3
# mpl.rcParams['boxplot.whiskerprops.color'] = 'C0'
mpl.rcParams['boxplot.capprops.linewidth'] = 3
mpl.rcParams['boxplot.medianprops.linewidth'] = 3
mpl.rcParams['boxplot.medianprops.color'] = bgcolor
# mpl.rcParams['boxplot.capprops.color'] = 'C0'
# mpl.rcParams['patch.facecolor'] = 'C0'


Path('../graficos').mkdir(exist_ok=True)

file_queue = [
#              'quantify_heads',
              'plot_distance_vs_corr',
              'plot_distance_vs_transcr',
#              'plot_motherlength_hist',
#              'plot_corr_motherlength',
#              'plot_transcription_hist',
#              'plot_sum_rc_by_head_bp',
              'plot_downstream_upstream',
              'plot_downstream_upstream_transcr',
              ]

if __name__ == '__main__':
    print('='*80, f"\n{'INICIANDO PLOTS':^80}\n",
          '='*80, '\n', '-'*80, '\n', sep='')

    for name in file_queue:
        print('='*80, f"\n{f'RUNNING {name}':^80}\n", '='*80, sep='')
        __file__ == name + '.py'
        exec('import ' + name)
        if show_flag:
            plt.show()
        plt.close('all')
