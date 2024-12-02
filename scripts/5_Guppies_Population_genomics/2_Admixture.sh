#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=64gb
#SBATCH -N 1
#SBATCH -n 20

# Input data
admixture=admixture_linux-1.3.0/admixture
input_bed=Guppies_prunned
# Number of runs as we are doing 10 runs per each K in total
run=$1

# Lowest K
Klow=1
# Maximum K
Khigh=15

# We create directory for given run
mkdir ${run}
cd ${run}
# We copy plink files to run directory
cp ${input_bed}.* ./
new_bed=$(basename ${input_bed})

for ((K=$Klow;K<=$Khigh;K++)); \
do
        ${admixture} -s ${RANDOM} -j19 --cv ${new_bed}.bed ${K} | tee log.${new_bed}.${K}.run${run}.out
        mv ${new_bed}.${K}.Q ${new_bed}.${K}.run${run}.Q
done
