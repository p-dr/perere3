# description: Removes non-gene annotations from S. mansoni genome annotations.
# in: pardir/'genome_annotation/schistosoma_mansoni.PRJEA36577.WBPS12.annotations.gff3'
# out: pardir/'genome_annotation/gene_annotations.gff3'

from pandas import read_csv
from utils import pardir, redo_flag
from wget import download

raw_annotations_path = pardir/'genome_annotation/schistosoma_mansoni.PRJEA36577.WBPS12.annotations.gff3'
annotations_path = pardir/'genome_annotation/gene_annotations.gff3'

if not raw_annotations_path.exists() or redo_flag:
    annotations_url = 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS13/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS13.annotations.gff3.gz'
    download(annotations_url, raw_annotations_path)
    

#============== REMOVER ANOTAÇÕES NÃO-GÊNICAS ===============#

print(f"Arquivos abertos. Lendo '{str(raw_annotations_path)}'...")
raw_annotations = read_csv(str(raw_annotations_path), sep='\t', skiprows=3, header=None, converters = {8: lambda s: 'gene_id='+s[s.find(':')+1:s.find(';')]})

print('Leitura encerrada. Removendo anotações não-gênicas...')

filtered_annotations = raw_annotations.loc[raw_annotations[2] == 'gene']
filtered_annotations[[3, 4]] = filtered_annotations[[3, 4]].astype(int)

filtered_annotations.to_csv(str(annotations_path), sep='\t', index=False, header=None)

print(f"Anotações gênicas mantidas em '{str(annotations_path)}'.")
