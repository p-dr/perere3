from pandas import read_csv
from utils import pardir, redo_flag
from matplotlib import pyplot

indir = pardir/'alinhamentos/SRA_vs_heads/only_mapped'
outpath = pardir/'genome_annotation/real_heads'

for sam_file in indir.glob('*'):
    
    
