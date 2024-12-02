#!/bin/bash

# Quality-based Missingness
# Here we change normal genotypes to missing if they fail the quality tresholds

#SBATCH --time=12:00:00
#SBATCH --mem=16gb
#SBATCH -o %j_Qual_Miss.out

module load python/3.9.2

# First input is the original file name and directory.
input_file=$1
echo '==========='
echo 'Input file:'
echo "$input_file"
# Second input is the final file name and directory.
output_file=$2
echo '==========='
echo 'Output file:'
echo "$output_file"

# Third input is a failed genotypes file which shows the failure reason.
fd_file=$3
echo '==========='
echo 'Failed genotypes file:'
echo "$fd_file"

# Fourth input is a file containing DP threshold for each sample (just two columns, min and max) in CORRECT ORDER (same as seen in VCF).
dp_file=$4
echo '==========='
echo 'Individuals DP file:'
echo "$dp_file"

# Fifth input is a GQ threshold (for example - 30).
echo '==========='
echo 'GQ threshold:'
echo "$GQ"
GQ=$5

# Sixth input is total lines count (for measuring progress sake).
total_lines=$(wc -l < ${input_file})

echo 'Starting at:'
date

python 4_Quality_based_missing.py \
--vcf $input_file \
--out $output_file \
--fd $fd_file \
--dp $dp_file \
--gq $GQ \
--lines ${total_lines}

echo 'Finished at:'
date
