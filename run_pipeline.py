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

u.logfile = u.LOG_PATH.open('a')
graph = nx.drawing.nx_pydot.read_dot(u.pardir/'pipeline')
start = "../genome_annotation/\n/all_together_now.tsv"

# argparser = argparse.ArgumentParser(
#     description='Run necessary scripts for creating ' + start.replace('/\n', ''))
# 
# argparser.add_argument('-f', '--startfrom',
#     default=None,
#     help='Start running pipeline from <script>')
# argparser.add_argument('-l', '--list', '--print',
#     action='store_true',
#     help='Only prints the pipeline to be run')
# argparser.add_argument('-e', '--exclude',
#     nargs='*',
#     help='Do not run selected scripts')


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


def run_script(name):
    if name == 'aggregate_data.py':  # find better way
        return
    u.print_header(f'Running {name}...', log=True)
    script_path = u.scripts_dir/name
    suffix = script_path.suffix

    if suffix == '.py':
        sname = name.strip('.py')
        try:
            exec(f'import {sname}; main_ret = {sname}.main()', globals())  # del after?

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
            u.logfile.write(line)
            print(line, end='')

        print('Child process returned', process.returncode)
        if process.returncode:
            exit()


def plot_all(confirm=False):
    if not confirm or u.user_confirms():
        for script in u.plot_scripts_dir.glob('*.py'):
            u.print_header(script.stem)
            sp.run(('python', str(script)))


def main(start=start):
    ps = parent_scripts(graph, start)
    # args = argparser.parse_args()

    # if args.startfrom is not None:
    #     ps = ps[ps.index(args.startfrom):]

    u.log('Iniciando sessão com os seguintes scripts:\n', *ps)
    finished = False  # finished downloading all accs.

    while not finished:
        for script in ps:
            # if args.list:
            #     print(script)
            # else:
            run_script(script)

            if main_ret == FINISHED_FLAG:
                finished = True

    u.print_header('agregando resultados', log=True)
    run_script('aggregate_data.py')

    if u.args.plot:
        plot_all()


if __name__ == '__main__':
    main()

u.logfile.close()
