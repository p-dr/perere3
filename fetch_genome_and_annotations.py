# description: Downloads genome and annotations file from "https://parasite.wormbase.org/Schistosoma_mansoni_prjea36577a" if files not found on disk.
# in: (cloud)
# out: pardir/'seqs/sm_genome.fa'
# out: pardir/'genome_annotation/sm_annotations.gff3'
from wget import download
from pathlib import Path
from subprocess import call
from utils import redo_flag

genome_out_dir = Path('../seqs')
genome_out_dir.mkdir(exist_ok=True)

annot_out_dir = Path('../genome_annotation')
annot_out_dir.mkdir(exist_ok=True)

genome_url = 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz'
annotations_url = 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3.gz'


for url, out_name, out_dir in ((genome_url, 'sm_genome.fa', genome_out_dir),
                               (annotations_url, 'sm_annotations.gff3', annot_out_dir)):

    out_path = out_dir/out_name
    raw_out_name = url.split('/')[-1].strip('.gz')
    raw_out_path = out_dir/raw_out_name

    if not out_path.exists() or redo_flag:
        print(f'Baixando {str(out_name)}')
        download(url, str(out_dir))
        call(f'gunzip {str(raw_out_path)}', shell=True)
        raw_out_path.rename(out_path)
    else:
        print(f'{out_path} já se encontra no disco.')

