#!/bin/bash

module load gatk/4.1.4.1

# Here we filter VCF files for correct SNPs.

#Population
p=$1

gatk --java-options "-Xmx16g" VariantFiltration \
-V 1_INFO_filtered_SNPs/INFO_unfiltered/${p}.vcf.gz \
-O 1_INFO_filtered_SNPs/${p}.vcf.gz \
--create-output-variant-index true \
--filter-expression "QD<2.0" --filter-name "QualityByDepth" \
--filter-expression "FS>60.0" --filter-name "FisherScore" \
--filter-expression "MQ<40.0" --filter-name "MappingQuality" \
--filter-expression "MQRankSum<-12.5" --filter-name "MappingQualityRankSum" \
--filter-expression "ReadPosRankSum<-8.0" --filter-name "ReadPositionRankSum"
