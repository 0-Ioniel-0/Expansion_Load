#!/bin/bash

################# BAM Manipulation #################
# BAM file produced in the previous step is still imperfect. It's best to remove not needed duplicates,
# which in result, makes files smaller. It's obligatory to add RG (read group) so that each sample has unique
# description and known sequencing run.

# Samples file
SAMPLE=$1
echo '========================'
echo 'Sample:'
echo "$SAMPLE"

# Providing directories #
IN=dir/to/${sample}.bam
# Picard tool
PIC=dir/to/tools/picard/picard.jar
TD=dir/to/temporary/files
# Library kit
library='unknown' # that's our example
# Sequencing platform
platform='unknown' # that's our example
unit='unit_1' # that's our example

echo '========================'
echo 'BAM refinement starts'
date
echo '===='
echo 'Name based sorting'
samtools sort -@ 5 -O bam -n -m 100M -o ${IN}.qsorted.bam ${IN}
# Note that this file is not indexed, as it CANNOT be done for query sorted files.
echo '===='
echo 'Fixing mates'
samtools fixmate -@ 5 -O bam ${IN}.qsorted.bam ${IN}.mates_fixed.bam 
echo '===='
echo 'Cooridnates based sorting'
samtools sort -@ 5 -m 100M -O bam -o ${IN}.csorted.bam ${IN}.mates_fixed.bam
samtools index ${IN}.csorted.bam
echo '===='
echo 'Marking, removing duplicates and adding Read Groups'
java -Xmx60G -jar ${PIC} MarkDuplicates \
        I=${IN}.csorted.bam \
        O=${IN}.rm_dup.bam \
        M=${input_bam}.txt \
        TMP_DIR=${TD} \
        CREATE_INDEX=true \
        REMOVE_DUPLICATES=true
java -Xmx60G -jar ${PIC} AddOrReplaceReadGroups \
        I=${IN}.rm_dup.bam \
        O=${DIR}/${sample}_ready.bam \
        RGLB=${library} \
        RGPL=${platform} \
        RGPU=${unit} \
        RGSM=${sample} \
        CREATE_INDEX=true        

# Remove intermediate files #
rm ${IN}.qsorted.bam
rm ${IN}.mates_fixed.bam
rm ${IN}.csorted.bam
rm ${IN}.csorted.bam.bai
rm ${IN}.rm_dup.bam
rm ${IN}.rm_dup.bam.bai

echo 'ALL finished'
date
