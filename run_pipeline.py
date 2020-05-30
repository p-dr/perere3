# pipeline = (
# gen_heads_complete_tr.py
# quantify_heads.py
# get_gene_annotations.py
# count_SRA_reads.py
# relate_heads_to_near_genes.py
#     correlate_heads_to_near_genes.py
#     aggregate_data.py
# )

import networkx as nx
import subprocess as sp
import utils as u
from fetch_SRA_data import FINISHED_FLAG
from sys import argv

graph = nx.drawing.nx_pydot.read_dot(u.pardir/'pipeline')
start = "../genome_annotation/\n/all_together_now.tsv"
# Which scripts must only be run after all pipeline loops.
run_after = ['correlate_heads_to_near_genes.py', 'aggregate_data.py']


def above_tree(graph, start):
    ret = []
    ret.append(list(graph.predecessors(start)))

    while True:
        new_line = []
        for parent_node in ret[-1]:
            new_line += list(graph.predecessors(parent_node))

        if not new_line:
            break
        ret.append(new_line)

    # parent_scripts = [[s for s in line if s[-3:] in {'.py', '.sh'}] for line in ret]
    # return [i for i in parent_scripts if i]
    return ret


def stratified_parent_scripts(graph, start, flag=True):
    # flag controls if scripts or outfiles.
    # False means types of start and outputs are
    # different, True means the opposite.
    ret = []
    old_line = set(graph.predecessors(start))

    while old_line:
        if flag and old_line:
            ret.append(old_line)

        new_line = set()
        for parent_node in old_line:
            new_line |= set(graph.predecessors(parent_node))

        old_line = new_line
        flag = not flag

    return ret


def parent_scripts(graph, start, flag=True):
    """Return possible scripts pipeline to generate start file.
       - flag controls if scripts or its outfiles are to be returned:
             True means types of start and outputs are
             different, False means the opposite.

      Future: parallelize scripts which can be run in parallel.
    """
    ret = []
    old_line = set(graph.predecessors(start))

    while old_line:
        if flag and old_line:
            ret = [i for i in ret if i not in old_line]
            ret = list(old_line) + ret

        new_line = set()
        for parent_node in old_line:
            new_line |= set(graph.predecessors(parent_node))

        old_line = new_line
        flag = not flag

    return ret


def run_script(name, finished=False):
    if (name in run_after) and (not finished):
        return
    u.print_header(f'Running {name}...', log=True)
    script_path = u.scripts_dir/name
    suffix = script_path.suffix

    if suffix == '.py':
        sname = name.strip('.py')
        try:
            exec(f'import {sname}; main_ret = {sname}.main(); del {sname}', globals())  # del after?

        except FileExistsError:
            u.log('FILEEXISTS')
            

    elif suffix == '.sh':
        process = sp.Popen(
            ['sh', name],
            stdout=sp.PIPE,
            stderr=sp.STDOUT,
            universal_newlines=True,
        )

        while process.poll() is None:
            line = process.stdout.readline()
            u.log(line.strip())

        if process.returncode:
            raise sp.CalledProcessError(f'Child process returned {process.returncode}.')


def plot_all(confirm=False):
    if not confirm or u.user_confirms():
        for script in u.plot_scripts_dir.glob('*.py'):
            u.print_header(script.stem)
            sp.run(('python', str(script)))


def run_end_scripts(scripts, plot=False):
    u.print_header('running end scripts...', log=True)
    u.log('End scripts:', *scripts, end='\t\n')

    for script in scripts:
        run_script(script, finished)

    if plot:
        u.print_header('plotting...', log=True)
        plot_all()


def main(start=start):

    local_plot_flag = u.args.plot
    u.args.plot = False
    u.plot_flag = False
    ps = parent_scripts(graph, start)

    u.log('Iniciando sessão com os seguintes script:\n', *ps)
    # Finished downloading all accs?
    finished = False

    while not finished:
        for script in ps:
            run_script(script)

            if main_ret == FINISHED_FLAG:
                finished = True  # save in this flag because main_ret will keep changing after.
    run_end_scripts(run_after, plot=local_plot_flag)

    # reset utils flags
    u.args.plot, u.plot_flag = 2 * [local_plot_flag]


if __name__ == '__main__':
    main()
