# description: Builds heads ht2 database.
# in: pardir/'seqs/heads.fa'
# out: pardir/'heads_ht2'
if [ ! -e ../heads_ht2 ]; then
    mkdir ../heads_ht2
    hisat2-build ../seqs/heads.fa ../heads_ht2/heads
else
    echo "heads' ht2 already exists."
fi
