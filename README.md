# Expansion Load
Scripts related to wild guppies expansion load estimation.

This repository contains scripts regarding:

    reads trimming
    reads alignment
    bam files manipulation
    variant calling
    filtering by quality, ancestral alleles and conservation score
    calculating genetic load
    estimating effective populations size and more...
    
!!! For details on the workflow of this project, please look into flowchart.png figure. !!!

Note that most of the scripts need to be edited to change paths to files and programs. DP thresholds, number of individuals, individuals' positions in the VCF files etc. should also be modified accordingly.
These scripts were run with:

    Python (3.9)
    FastQC (0.11.8)
    Trimmomatic (0.39)
    bwa (0.7.10)
    samtools (1.6.0)
    bcftools (1.9)
    gatk (4.1.4.1)
    vcftools (0.1.14)
    R(4.3.3)
    pixy(1.2.7.beta1)
    admixture(1.3.0)
    SNPeff(5.0)
    plink(1.90)

For help or comments, contact: k.burda[a]amu.edu.pl
