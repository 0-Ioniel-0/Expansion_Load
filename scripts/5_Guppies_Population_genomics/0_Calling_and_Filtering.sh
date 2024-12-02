#!/bin/bash
#SBACTH -N 1
#SBATCH -n 5
#SBATCH --time=24:00:00
#SBATCH --mem=32G

module load samtools/1.20
module load bcftools/1.20

# BAM files are already prepared in 1_Guppies_Mapping_and_Calling
# Joining many BAMs into one BCF
samtools mpileup -Ou -f Guppy.fasta -a AD -b BamFiles.txt -o Raw.Guppies.bcf

# Calling SNPs
bcftools call -mv -Ob -o Guppy.SNPS.bcf Raw.Guppies.bcf
bcftools index Guppy.SNPS.bcf

# Filtering SNPs
bcftools view -Oz -v snps -i 'MQ>40 && INFO/DP<5200 && INFO/DP>1300' -e 'QUAL <= 30 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || SCBZ > 6' -o Guppies_filtered.vcf.gz Guppy.SNPS.bcf
bcftools index Guppies_filtered.vcf.gz
