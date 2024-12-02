#!/bin/bash

# Calculating nucleoide diversity with pixy.
# Population
p=$1

pixy --stats pi \
--vcf ${p}.vcf.gz \
--populations pop/${p}.pop \
--window_size 1000000 \
--n_cores 10 \
--output_folder pixy_out \
--output_prefix ${p}
