from run_pipeline import run_script
import gen_heads
import utils as u
from count_SRA_reads import out_dir as count_dir

pipeline = (
    'gen_heads.py',
    'quantify_heads.py',
    'gen_complement_annotations.py',
#    'relate_heads_to_near_genes.py',
#    'build_heads_ht2.sh',
)


def main():
    old_redo = u.args.redo
    u.args.redo = True
    u.redo_flag = True

    u.clean(*count_dir.glob('*head*'), force=True)
    for script in pipeline:
        run_script(script)

    u.args.redo = old_redo
    u.redo_flag = old_redo


if __name__ == '__main__':
    main()
