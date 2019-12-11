from wget import download
from pathlib import Path

out_dir = Path('../seqs')
out_dir.mkdir(exist_ok=True)

genome_url = 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz'

annotations_url = 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3.gz'


for url in (genome_url, annotations_url):
	download(url, out_dir)
