# description: Builds S. mansoni genome ht2 database.
# in: pardir/'seqs/sm_genome.fa'
# out: pardir/'genome_ht2'

# Must be run from scripts dir.
if [ ! -e '../genome_ht2' ]; then
    mkdir ../genome_ht2
    hisat2-build ../seqs/sm_genome.fa ../genome_ht2/sm_genome
else
    echo 'Genome ht2 already exists.'
fi
