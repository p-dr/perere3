from pandas import read_csv
from utils import pardir

raw_annotations_path = pardir/'genome_annotation'/'schistosoma_mansoni.PRJEA36577.WBPS12.annotations.gff3'
annotations_path = (pardir/'genome_annotation'/'gene_annotations.gff3')


#============== REMOVER ANOTAÇÕES NÃO-GÊNICAS ===============#

print(f"Arquivos abertos. Lendo '{str(raw_annotations_path)}'...")
raw_annotations = read_csv(str(raw_annotations_path), sep='\t', skiprows=3, header=None, converters = {8: lambda s: 'gene_id='+s[s.find(':')+1:s.find(';')]})

print('Leitura encerrada. Removendo anotações não-gênicas...')

filtered_annotations = raw_annotations.loc[raw_annotations[2] == 'gene']
filtered_annotations[[3, 4]] = filtered_annotations[[3, 4]].astype(int)

filtered_annotations.to_csv(str(annotations_path), sep='\t', index=False, header=None)

print(f"Anotações gênicas mantidas em '{str(annotations_path)}'.")
