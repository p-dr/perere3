from utils import pardir
from os.path import getctime
from graphviz import Digraph
from matplotlib.pyplot import get_cmap
from matplotlib.colors import to_hex
from pprint import pprint

CMAP = get_cmap('cool').reversed()
pipeline = Digraph()
scripts_dir = pardir/'scripts'
ignore_list = ['*.sw*']
ext = ['.py', '.sh']


def path2str(p):
    if isinstance(p, str):
        return p
    # WRN: .replace of pathlib moves the file!
    return str(p.resolve()).replace(str(pardir), '..').replace('/','/\n/').replace('../\n', '..')


def eval_path(path_var, script_path):
    if path_var.startswith('('):
        return path_var
    elif script_path.suffix != '.py':
        return eval(path_var)
    exec(f'import {script_path.stem} as s; eval_path.ret = s.{path_var}')
    return eval_path.ret


def parse_script(script_path):
    script_tags = {}

    with script_path.open('r') as script_file:
        print('Processing', script_path)

        for line in script_file:
            if line.startswith('#'):

                tag = line.strip('#').strip().split(': ')

                if len(tag) == 1:
                    continue
                else:
                    tag, value = tag
                    tag = tag.strip()

                if tag == 'description':
                    script_tags[tag] = value.strip()

                elif tag in ('in', 'out'):
                    value = value.split()
                    value = [eval_path(path, script_path) for path in value]
                    script_tags[tag] = script_tags.get(tag, []) + value

            elif line != '\n':
                break

    return script_tags


def build_edges(script_tags, script_name, file2color_map):
    pipeline.attr('node', color='black', shape='box', style='solid', tooltip=script_tags['description'])
    pipeline.node(script_name)

    pipeline.attr('node', shape='folder', tooltip='File or folder.', style='filled')

    for inpath in script_tags.get('in', ()):
        pipeline.attr('node', color=file2color_map[inpath])
        pipeline.edge(path2str(inpath), script_name)

    for outpath in script_tags.get('out', ()):
        pipeline.attr('node', color=file2color_map[outpath])
        pipeline.edge(script_name, path2str(outpath))


def ctimes2color(all_iofiles_ctime):
    """Returns dictionary with colors from dictionary of creation times in seconds."""
    ALPHA = '66'
    ret = {'(cloud)': 'white'}
    ctimes = sorted(all_iofiles_ctime.values())
    color_step = CMAP.N // (len(ctimes)-1)  # why must it be integer?

    for io_str_path, ctime in all_iofiles_ctime.items():
        color = to_hex(CMAP(ctimes.index(ctime) * color_step)) + ALPHA
        ret.update({io_str_path: color})

    return ret


def main():

    untagged_scripts = []
    all_iofiles_ctime = {}
    all_scripts_tags = {}

    for script_path in (p for p in scripts_dir.glob('*') if p.suffix in ext):

        script_name = script_path.name
        script_tags = parse_script(script_path)

        if script_tags:
            all_scripts_tags.update({script_name: script_tags})

            for io_path in script_tags.get('in', []) + script_tags.get('out', []):

                if isinstance(io_path, str) or any(io_path.match(patt) for patt in ignore_list):
                    continue

                elif not io_path.exists():
                    # Future: annotate it does not exist.
                    creation_time = 0

                elif io_path.is_dir():

                    dir_children = list(io_path.glob('*'))

                    if dir_children:
                        creation_time = max(map(getctime, dir_children))
                    else:
                        creation_time = 0

                else:
                    creation_time = io_path.stat().st_mtime

                all_iofiles_ctime[io_path] = all_iofiles_ctime.get(io_path, creation_time)

        else:
            untagged_scripts.append(script_name)

    file2color_map = ctimes2color(all_iofiles_ctime)

    for script_name, script_tags in all_scripts_tags.items():
        build_edges(script_tags, script_name, file2color_map)


    print('Não utilizados:', *untagged_scripts)
    pipeline.render(pardir/'pipeline', format='svg', cleanup=True)
    pipeline.render(pardir/'pipeline', format='pdf')#, cleanup=True)#, view=True)


if __name__ == '__main__':
    main()
