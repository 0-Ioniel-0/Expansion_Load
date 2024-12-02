import argparse
import os
from collections import namedtuple
from typing import Set

# Define a named tuple for BED entry
BedEntry = namedtuple("BedEntry", ['chrom', 'start'])


def read_bed_file(bed_file):
    bed_data: set[tuple] = set()
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start = line.strip().split('\t')[:2]
            bed_data.add(BedEntry(chrom=chrom, start=start))
    return bed_data


def filter_vcf_with_bed_info(vcf_file, bed_data, output_file):
    with open(vcf_file, 'r') as vcf, open(output_file, 'w') as out:
        for line in vcf:
            if line.startswith('#'):
                out.write(line)
            else:
                fields = line.strip().split('\t')
                chrom, pos = fields[0], fields[1]

                if BedEntry(chrom=chrom, start=pos) in bed_data:
                    out.write(line)
                else:
                    print(f"Variant not found in BED: {chrom}:{pos}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FILTERING 4fds - Here we filter VCF file by 4-fold degenerated "
                                                 "sites.")
    parser.add_argument('-v', '--vcf', help='Input vcf file. Must be uncompressed.', type=str, required=True)
    parser.add_argument('-o', '--out', help='Output vcf file.', type=str, required=True)
    parser.add_argument('-b', '--bed', help='BED file with 4fds', type=str, required=True)
    args = parser.parse_args()

    print("Hello!")
    print("We are running FILTERING 4-fold-degenerated-sites script")
    print("Input File: ", args.vcf)
    print("Output File: ", args.out)
    print("Working...")

    if not (os.path.exists(args.vcf) and os.path.exists(args.bed)):
        print("Error: One or more input files do not exist.")
    else:
        BED_data = read_bed_file(args.bed)
        filter_vcf_with_bed_info(args.vcf, BED_data, args.out)

    print("Finished.")

