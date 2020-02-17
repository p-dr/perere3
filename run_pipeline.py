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
from utils import pardir, scripts_dir, LOG_PATH
from sys import argv

log_file = LOG_PATH.open('a')
graph = nx.drawing.nx_pydot.read_dot(pardir/'pipeline')
start = "pardir/\n/'genome_annotation/\n/all_together_now.tsv'"

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
    script_path = scripts_dir/name
    suffix = script_path.suffix
    interpreter = {'.py': 'python', '.sh': 'sh'}
    command = (interpreter[suffix], str(script_path))
    print('Running', name, '...', end=' ')
    sp.run(command, stdout=log_file, stderr=log_file)
    print('Done.')


def main(start=start):
    ps = parent_scripts(graph, start)
   
    if '--from' in argv:
        start_script = argv[argv.index('--from') + 1]
        ps = ps[ps.index(start_script):]

    for script in ps:
        if '--list' in argv:
            print(script)
        else:
            run_script(script)


if __name__ == '__main__':
    main()

log_file.close()
