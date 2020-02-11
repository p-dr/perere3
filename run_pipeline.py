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
from utils import pardir
from sys import argv

graph = nx.drawing.nx_pydot.read_dot(pardir/'pipeline')
start = 'quantify_heads.py'

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


def parent_scripts(graph, start, flag=True):
    # flag controls if scripts or outfiles.
    ret = []
    old_line = list(graph.predecessors(start))

    while old_line:
        new_line = []
        for parent_node in old_line:
            new_line += list(graph.predecessors(parent_node))

        if flag and new_line:
            ret.append(new_line)

        old_line = new_line
        flag = not flag

    return ret

print(parent_scripts(graph, "pardir/\n/'genome_annotation/\n/all_together_now.tsv'"))
