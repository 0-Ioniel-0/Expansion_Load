#!/bin/bash

# GVCF to VCF
# Having some poor variants removed, let's genotype all remaining sites.

#SBATCH --time=3:00:00
#SBATCH --mem=16gb
#SBATCH -o %j_gvcf-vcf.out

module load gatk/4.1.4.1

# First input is the original file name and directory.
input_file=$1
echo 'Input file:'
echo '$input_file'
# Second input is the final file name and directory.
output_file=$2
echo 'Output file:'
echo '$output_file'

echo 'Starting'
date

gatk --java-options "-Xmx16G" GenotypeGVCFs \
-R reference.fasta \
-V $input_file \
-O $output_file \
-all-sites

date
echo 'Finished'