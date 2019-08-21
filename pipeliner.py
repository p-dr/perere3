from utils import pardir
from os.path import getctime
from graphviz import Digraph
from matplotlib.pyplot import get_cmap
from matplotlib.colors import to_hex
from pprint import pprint

CMAP = get_cmap('cool').reversed()
pipeline = Digraph()
scripts_dir = pardir/'scripts'


def parse_script(script_path):

    script_tags = {}

    with script_path.open('r') as script_file:

        for line in script_file:
            if line.startswith('#'):

                tag = line.strip('#').strip().split(': ')

                if len(tag) == 1:
                    continue

                if tag[0] != 'description':
                    tag[1] = tag[1].split(' ')

                script_tags.update([tag])
                
            elif line != '\n':
                break

    return script_tags


def build_edges(script_tags, script_name):
    pipeline.attr('node', color='black', shape='box', style='solid', tooltip=script_tags['description'])
    pipeline.node(script_name)

    pipeline.attr('node', shape='folder', tooltip='File or folder.', style='filled')
    
    for inpath in script_tags.get('in', ()):
        pipeline.attr('node', color=file2color_map[inpath])
        pipeline.edge(inpath.replace('/','/\n/'), script_name)

    for outpath in script_tags.get('out', ()):
        pipeline.attr('node', color=file2color_map[outpath])
        pipeline.edge(script_name, outpath.replace('/','/\n/'))


def ctimes2color(all_iofiles_ctime):
    """Returns dictionary with colors from dictionary of creation times in seconds."""
    ALPHA = '66'
    ret = {'(cloud)': 'white'}
    ctimes = sorted(all_iofiles_ctime.values())
    color_step = CMAP.N // (len(ctimes)-1)  # why must it be integer?
    for io_str_path, ctime in all_iofiles_ctime.items():
        # trim out 'pardir/"' and '"'
        io_path = eval(io_str_path)
        color = to_hex(CMAP(ctimes.index(ctime) * color_step)) + ALPHA
        ret.update({io_str_path: color})

    return ret


if __name__ == '__main__':

    untagged_scripts = []
    all_iofiles_ctime = {}
    all_scripts_tags = {}
    
    for script_path in scripts_dir.glob('*[.py,.sh]'):

        script_name = script_path.name
        script_tags = parse_script(script_path)

        if script_tags:
            all_scripts_tags.update({script_name: script_tags})

            for io_str_path in script_tags.get('in', []) + script_tags.get('out', []):
                if io_str_path != '(cloud)':
                    # trim out 'pardir/"' and '"'
                    io_path = pardir/io_str_path[8:-1]

                    if io_path.is_dir():

                        dir_children = list(io_path.glob('*'))

                        if dir_children:
                            creation_time = max(map(getctime, dir_children))
                        else:
                            creation_time = 0
                            
                    else:
                        creation_time = io_path.stat().st_mtime

                    all_iofiles_ctime[io_str_path] = all_iofiles_ctime.get(io_str_path, creation_time)
                
        else:
            untagged_scripts.append(script_name)

    file2color_map = ctimes2color(all_iofiles_ctime)
    
    for script_name, script_tags in all_scripts_tags.items():
        build_edges(script_tags, script_name)
        

    print('Não utilizados:', *untagged_scripts)
    pipeline.render(pardir/'pipeline', format='svg', cleanup=True)
    pipeline.render(pardir/'pipeline', format='pdf', cleanup=True)#, view=True)
