from matplotlib.pyplot import close
from importlib import reload
from pathlib import Path

Path('../graficos').mkdir(exist_ok=True)

file_queue = {'plot_distance_vs_corr',
              'plot_distance_vs_transcr',
              'plot_corr_motherlength',
              'plot_transcription_hist',
              'plot_sum_rc_by_head_bp',
             }

if __name__ == '__main__':
    for name in file_queue:
        print('='*80, f"\n{f'RUNNING {name}':^80}\n", '='*80, sep='')
        __file__ == name + '.py'
        exec('import ' + name)
        close()

