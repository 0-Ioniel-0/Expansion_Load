#!/bin/bash

#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time=01:00:00
#SBATCH -o slurm_out/%j_pre-filtering.out

module load bcftools/1.9

#Pre-filtering with max 30% missing data
#This step is important as following ones are with custom scripts so the less poor data, the better.
#Chromosome and population
c=$1
p=$2

echo $c
echo $p
echo 'Starting'
date

bcftools view --threads 5 -i 'F_MISSING <= 0.3' \
-O z \
-o 1_30_missing_pre_filtering/${c}_${p}.g.vcf.gz \
~/${c}_${p}.g.vcf.gz

date
echo 'Finished'
