#!/bin/bash

#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time=01:00:00
#SBATCH -o slurm_out/%j_pre-filtering.out

module load bcftools/1.9

#Pre-filtering with max 30% missing data based on minimal DP and GQ (treated as missing if not meeting thresholds).
#This step is important as following ones are with custom scripts so the less poor data, the better.

input=$1
output=filtered_$(basename ${input} .g.vcf.gz).vcf
echo 'Starting'
date

bcftools view --threads 5 -i '((F_PASS(FMT/GQ>=30 & FMT/DP>=6 & GT!="mis") > 0.7) & (TYPE="snp")) | ((F_PASS(FMT/RGQ>=30 & FMT/DP>=6 & GT!="mis") > 0.7) & (TYPE="ref"))' -Ov -o ${output} ${input}

date
echo 'Finished'
