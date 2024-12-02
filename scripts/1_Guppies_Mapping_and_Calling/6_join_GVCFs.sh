#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --mem=64GB

################# GATK gVCF combining #################
# gVCFs produced in previous step should be now combined to allow for subsequent joined calling.
# In this example, we have a population with 14 samples.

# Population name, given first
POPULATION=$1
echo '========================'
echo 'Population:'
echo "$POPULATION"
# Chromosome, given second
CHROMOSOME=$2
echo '========================'
echo 'Chromosome:'
echo "$POPULATION"
# Samples file (of given population [with the same name], horizontal list of names, in our example - 14)
samples=${POPULATION}.txt
s1=$(echo "$samples" | cut -f1)
s2=$(echo "$samples" | cut -f2)
s3=$(echo "$samples" | cut -f3)
s4=$(echo "$samples" | cut -f4)
s5=$(echo "$samples" | cut -f5)
s6=$(echo "$samples" | cut -f6)
s7=$(echo "$samples" | cut -f7)
s8=$(echo "$samples" | cut -f8)
s9=$(echo "$samples" | cut -f9)
s10=$(echo "$samples" | cut -f10)
s11=$(echo "$samples" | cut -f11)
s12=$(echo "$samples" | cut -f12)
s13=$(echo "$samples" | cut -f13)
s14=$(echo "$samples" | cut -f14)

# We used GATK 4.1.4.1
GA=dir/to/gatk

# Reference genome file, given third
REF=$3

# Specify number of threads available (n), PE optiom idicates pair-end data
echo '========================'
echo 'GVCF joining starts'
date

$GATK --java-options "-Xmx64G" CombineGVCFs \
-R ${RE} \
-O ${CHROMOSOME}_${POPULATION}.g.vcf.gz \
--variant ${CHROMOSOME}_${s1}.g.vcf.gz \
--variant ${CHROMOSOME}_${s2}.g.vcf.gz \
--variant ${CHROMOSOME}_${s3}.g.vcf.gz \
--variant ${CHROMOSOME}_${s4}.g.vcf.gz \
--variant ${CHROMOSOME}_${s5}.g.vcf.gz \
--variant ${CHROMOSOME}_${s6}.g.vcf.gz \
--variant ${CHROMOSOME}_${s7}.g.vcf.gz \
--variant ${CHROMOSOME}_${s8}.g.vcf.gz \
--variant ${CHROMOSOME}_${s9}.g.vcf.gz \
--variant ${CHROMOSOME}_${s10}.g.vcf.gz \
--variant ${CHROMOSOME}_${s11}.g.vcf.gz \
--variant ${CHROMOSOME}_${s12}.g.vcf.gz \
--variant ${CHROMOSOME}_${s13}.g.vcf.gz \
--variant ${CHROMOSOME}_${s14}.g.vcf.gz

echo 'ALL finished'
date
