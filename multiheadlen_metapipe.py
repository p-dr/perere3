import run_pipeline
import gen_heads
import utils as u
import shutil
from count_SRA_reads import out_dir as count_dir
from aggregate_data import outpath as outagg

lengths = (50, 100, 150, 200, 300, 400, 500, 700, 1000, 1250, 1500, 2000)
dir_suffix = 'bp_counted_reads'
u.args.clean = True


def main():
    for length in lengths:
        out_dir = u.pardir/f'{length}{dir_suffix}'
        dirprefix = out_dir.stem.split('_')[0]

        if out_dir.exists():
            print(out_dir, 'already exists.')
        else:
            u.print_header(length, 'bp pipeline running', log=True)
            gen_heads.HEAD_LEN = length
            run_pipeline.run_script('remake_heads_pipeline.py')
            run_pipeline.main()

            u.log('Saving results...')
            shutil.copytree(count_dir, out_dir)
            shutil.copy(outagg, outagg.parent/(dirprefix + '_' + outagg.name))


if __name__ == '__main__':
    main()
