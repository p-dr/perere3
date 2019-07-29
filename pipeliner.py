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

    untagged_scripts = []
    
    for script_path in scripts_dir.glob('*[.py,.sh]'):
        script_name = script_path.name
        script_tags = parse_script(script_path)

        if not script_tags:
            untagged_scripts.append(script_name)
            continue

        pipeline.attr('node', color='lightgray', shape='box', style='filled', tooltip=script_tags['description'])
        pipeline.node(script_name)

        pipeline.attr('node', color='black', shape='folder', style='solid', tooltip='File or folder.')
        
        for inpath in script_tags.get('in', ()):
            pipeline.edge(inpath, script_name)

        for outpath in script_tags.get('out', ()):
            pipeline.edge(script_name, outpath)


    print('Não utilizados:', *untagged_scripts)
    pipeline.render(pardir/'pipeline', format='svg', cleanup=True)
    pipeline.render(pardir/'pipeline', format='pdf', cleanup=True)#, view=True)
