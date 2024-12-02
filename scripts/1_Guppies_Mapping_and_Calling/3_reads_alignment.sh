#!/bin/bash

################# ALIGNMENT ################# 
# Having trimmed reads and indexed reference genome, one can make a BAM file.

# Input reads
IN=dir/to/TRIMMED
# Out BAMs
OUT=dir/to/BAMs
# Reference files. Index files should be mentioned only with directory and prefix, without .fasta !!!
R=dir/to/ref
# Temporary directory
T=dir/to/temp
# Sample name given by the outside loop. Each sample was run separately in the same time.
SAMPLE=$1
echo '========================'
echo 'Sample:'
echo "$SAMPLE"

echo '========================'
echo 'Aligning starts'
date
# Aligning reads with bwa mem. To save space on the disc, it's best to forward output to samtools view which
# changes SAM (created by bwa) into BAM.
bwa mem ${R} ${IN}/${SAMPLE}_1P.fastq.gz ${IN}/${SAMPLE}_2P.fastq.gz |\
samtools view -b -F 4 - > ${OUT}/${SAMPLE}.bam
# ^ option -F 4 makes sure only mapped reads are saved.
echo '========================'
echo 'Aligning finished'
date
