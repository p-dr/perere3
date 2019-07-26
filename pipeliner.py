from utils import pardir
from graphviz import Digraph

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

                tag[1] = tag[1].replace('/','/\n/')
                
                if tag[0] != 'description':
                    tag[1] = tag[1].split(' ')

                script_tags.update([tag])
                
            elif line != '\n':
                break

    return script_tags


if __name__ == '__main__':
    for script_path in scripts_dir.glob('*[.py,.sh]'):
        script_name = script_path.name
        script_tags = parse_script(script_path)

        # if not script_tags: continue
        
        pipeline.attr('node', color='lightgray', shape='box', style='filled', tooltip='hooooooo')
        pipeline.node(script_name)

        pipeline.attr('node', color='black', shape='folder', style='solid')
        
        for inpath in script_tags.get('in', ()):
            pipeline.edge(inpath, script_name)

        for outpath in script_tags.get('out', ()):
            pipeline.edge(script_name, outpath)

    pipeline.render(pardir/'pipeline', format='pdf', cleanup=True, view=True)
