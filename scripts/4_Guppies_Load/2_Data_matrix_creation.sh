#!/bin/bash

#SBATCH --time=02:00:00

module load python/3.10

# Here, we check for load in all individuals within a population VCF file.
# The final matrix will be huge. It will have rows with each derived alleles, with columns: Sample Chromosome Region CS Site_Load Site_Effect(for LOFs) Allele_Count

# Population
p=$1

echo $1
echo $2

date

python 2_Data_matrix_creation.py --vcf 1_VCFs/${c}_${p}.vcf --out 2_Load --pop $p
date
