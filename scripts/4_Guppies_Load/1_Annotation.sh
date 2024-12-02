#!/bin/bash

module load java8

# First one needs to look for the right database
java -jar snpEff.jar databases | grep 'Guppy'

# Now, let's annotate. In our case it was important to consider scaffolds names (NCBI vs. GenBank ones).
java -jar snpEff.jar Guppy_female_1.0_MT.99 Guppy.vcf > Annotated_Guppy.vcf
