#!/bin/bash

#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time=24:00:00
#SBATCH --mem=64g

################# GATK variant calling #################
# Let's call all sites in GATK, a genomic VCF (gVCF) is firstly produced for each individual separately.
# To make things faster it's good to call chromosomes, or even parts of chromsomes separately (-L option).

# Sample file
SAMPLE=$1
echo '========================'
echo 'Sample:'
echo "$SAMPLE"
# Reference genome file
REF=$2
# Chromosome to call
CHROMOSOME=$3
# Output file
OUT=output_directory/${CHROMOSOME}_$(basename $SAMPLE .bam).g.vcf.gz

# We used GATK 4.1.4.1
GATK=dir/to/gatk

echo '========================'
echo 'Calling starts'
date

# REMINDER - we are calling ALL sites! They are going to be needed in neutral nucleotide diversity calculation.
$GATK --java-options "-Xmx64g" HaplotypeCaller \
-R $REF \
-I $SAMPLE \
-O $OUT \
-ERC BP_RESOLUTION \
--pcr-indel-model NONE \
--output-mode EMIT_ALL_CONFIDENT_SITES \
--native-pair-hmm-threads 20 \
-L $CHROMOSOME

echo 'Calling finished'
date
