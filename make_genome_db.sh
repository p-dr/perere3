# description: Makes BLAST genome database.
# in: pardir/'seqs/sm_genome.fa'
# out: pardir/'sm_genome_db'

# Must run from scripts directory.
if [ ! -e '../sm_genome_db' ]
    then
        makeblastdb -in ../seqs/sm_genome.fa -out ../sm_genome_db/sm_genome_db -dbtype nucl -title "Schistosoma mansoni's genome"
else
    echo 'Genome BLAST database already exists.'
fi

