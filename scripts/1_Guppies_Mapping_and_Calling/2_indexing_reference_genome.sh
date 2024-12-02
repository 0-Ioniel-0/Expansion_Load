#!/bin/bash

################# Indexing Reference Genome #################
# In order to map reads to the reference genome, 
# it must be indexed first with the same bwa, alignment is going to be performed with.

# Reference genome file
REF=/dir/to/reference/ref.fasta

echo '========================'
echo 'Indexing starts'
date
# Guppy genome is ~700MB so it's big enough to index it with a 'bwtsw' method designed for big genomes.
bwa index -p Prefix -a bwtsw ${REF}
# Variant calling step also requires index, but this time it's a different one, you can create it with samtools.
samtools faidx ${REF}
echo 'Indexing finished'
date
