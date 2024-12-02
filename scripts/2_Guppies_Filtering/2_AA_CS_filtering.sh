#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=16gb
#SBATCH -o slurm_out/%j_AA_CS.out

# AA and CS filtration
# Having some poor variants removed, let's add info on AA and CS.

module load python/3.9.2

# First input is the original file (uncompressed!) name and directory.
input_file=$1
echo '==========='
echo 'Input file:'
echo "${input_file}"
# Second input is the final file name and directory.
output_file=output_dir/$(basename ${input_file})
echo '==========='
echo 'Output file:'
echo "${output_file}"

# Third input is a BED file with ancestral alleles in 4th column and conservation score in 5th (after chromosome, start and end.
# IMPORTANT NOTICE - this script works with 1-based BED, so NOT normal BED, which is 0-based. You can change your BED to 1-based
# by adding 1 to both START and END. Also, it's best to generally run this script on one chromosome, using just single chromosome
# BED file as well.
bed_file=$2
echo '==========='
echo 'BED file:'
echo "${bed_file}"

# Fourth input is a name for failed sites file, wich shows the reason for failure, e.g. no AA or no CS.
fd_file=3_AA_CS_filtered_sites/parts/second_parts/$(basename ${input_file}).fd
echo '==========='
echo 'Failed sites file:'
echo "${fd_file}"

# Fifth input is total lines count (for measuring progress sake).
total_lines=$(wc -l < ${input_file})

echo 'Starting'
date

python 3_AA_CS_filtered_sites.py \
--vcf "${input_file}" \
--out "${output_file}" \
--bed "${bed_file}" \
--fd "${fd_file}" \
--lines "${total_lines}"

date
echo 'Finished'
