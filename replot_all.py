import matplotlib.pyplot as plt
import matplotlib as mpl
from importlib import reload
from pathlib import Path

plt.style.use('seaborn')
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#22848e',
                                                    '#f88904',
                                                    '#78c18c'])
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

mpl.rcParams['boxplot.patchartist'] = True

Path('../graficos').mkdir(exist_ok=True)

file_queue = ['plot_distance_vs_corr',
              'quantify_heads',
              'plot_motherlength_hist',
              'plot_distance_vs_transcr',
              'plot_corr_motherlength',
              'plot_transcription_hist',
              'plot_sum_rc_by_head_bp',
             ]

if __name__ == '__main__':
    for name in file_queue:
        print('='*80, f"\n{f'RUNNING {name}':^80}\n", '='*80, sep='')
        __file__ == name + '.py'
        exec('import ' + name)
        plt.close()

