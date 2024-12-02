#!/bin/bash

# FINAL MISSINGNESS FILTERING
# Let's filter sites with missingness fraction of max 30%.

#SBATCH --time=01:00:00
#SBATCH -o %j_filtering.out

module load vcftools
module load tabix

echo 'Starting at:'
date

# First input is the original file name and directory.
input_file=$1
echo '============'
echo 'Input file:'
echo "$input_file"
# Second input is the final file name and directory.
output_file=$2
echo '============'
echo 'Output file:'
echo "$output_file"

# Filtering
vcftools --max-missing 0.7 --vcf $input_file \
--recode --recode-INFO-all --out $output_file

# Indexing
tabix $output_file

echo 'Finished at:'
date
