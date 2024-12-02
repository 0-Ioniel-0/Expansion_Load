#!/bin/bash

module load plink/1.90

# Plink files manipulation

# Turning VCF into plink formats, prunning it and doing PCA
plink --vcf Guppies.vcf --keep-allele-order --make-bed --out Guppies
plink --bfile Guppies --indep-pairwise 200kb 1 0.5 --out prunned_snps
plink --bfile Guppies --extract pruned_snps.prune.in --make-bed --out Guppies_prunned
plink --bfile Guppies_pruned --pca --out Guppies_PCA

