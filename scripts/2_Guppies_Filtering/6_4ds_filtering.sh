#!/bin/bash

#SBATCH --time=0:15:00
#SBATCH --mem=8gb

#Filter genomic VCF for trully synonymous sites (4-fold degenerated sites).

date

input_file=$1
output_file=$(basename ${input_file})
bed_file=$2

module load python/3.9.2

python 6_4fds_filtering.py \
-v ${input_file} \
-o 6_4fds_filtering/${output_file} \
-b ${bed_file}

date
