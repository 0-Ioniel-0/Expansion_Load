import argparse

parser = argparse.ArgumentParser(description="Adding Ancestral Alleles and Conservation Scores - Here we filter VCF "
                                             "file by ancestral alleles and conservation score presence. Variants that "
                                             "have AA as either REF or ALT are outputted.")
parser.add_argument('-v', '--vcf', help='Input vcf file. Must be uncompressed.', type=str, required=True)
parser.add_argument('-o', '--out', help='Output vcf file.', type=str, required=True)
parser.add_argument('-b', '--bed', help='BED file with ancestral alleles in 4th column and conservation'
                                        'score in 5th (after chromosome, start and end.', type=str,
                    required=True)
parser.add_argument('-f', '--fd', help='Sites info file.', type=str, required=True)
parser.add_argument('-l', '--lines', help='Total lines in file (used for counting progress.', type=str,
                    required=True)
args = parser.parse_args()

print("Hello!")
print("We are running FILTERING Ancestral-Alleles/Conservation-Score script")
print("Input File: ", args.vcf)
print("Output File: ", args.out)
print("Sites info file: ", args.fd)
print("Working...")


def read_bed_file(bed_file):
    aa_bed_data = {}
    cs_bed_data = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end, allele, cs = line.strip().split('\t')
            aa_bed_data[(chrom, str(start))] = allele
            cs_bed_data[(chrom, str(start))] = cs
    return aa_bed_data, cs_bed_data


def update_vcf_with_bed_info(vcf_file, aa_bed_data, cs_bed_data, output_file, fd_file):
    with open(vcf_file, 'r') as vcf, open(output_file, 'w') as out:
        line_number = 1
        total_lines = int(args.lines)
        for line in vcf:

            reason = "None."
            chrom = "0"
            pos = "0"
            if line.startswith('#'):
                if line.startswith('##INFO=<ID=ReadPosRankSum'):
                    out.write(line)
                    out.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele.">' + '\n')
                    out.write('##INFO=<ID=CS,Number=1,Type=String,Description="Conservation score.">' + '\n')
                    out.write('##INFO=<ID=SwappedGenotype,Number=1,Type=String,Description="Information '
                              'on GT swap Ancestral/Derived when ALT is Ancestral.">' + '\n')
                else:
                    out.write(line)  # Write header lines directly to output file
            else:
                fields = line.strip().split('\t')
                chrom, pos, ref, alt = fields[0], str(fields[1]), fields[3], fields[4]
                info = fields[7]
                with open(fd_file, 'a') as fd:
                    fd.write(chrom + pos + "\t")
                print("Line", line_number, "of", total_lines, (line_number*100)/total_lines, "%")

                if (chrom, pos) in aa_bed_data:
                    allele_from_bed = aa_bed_data[(chrom, pos)]
                    cs_from_bed = cs_bed_data[(chrom, pos)]

                    if allele_from_bed in [ref, alt]:
                        samples = fields[9:]
                        if allele_from_bed == alt:
                            reason = "Swapped genotypes."
                            new_alt = str(ref)
                            new_ref = str(alt)
                            fields[3] = new_ref
                            fields[4] = new_alt
                            samples = fields[9:]
                            for i in range(len(samples)):
                                sample = samples[i]
                                GT = sample.split(':')[0]
                                if GT == "0/0" or GT == "0|0":
                                    GT = "1/1"
                                elif GT == "1/1" or GT == "1|1":
                                    GT = "0/0"
                                samples[i] = ':'.join([GT] + sample.split(':')[1:])
                        fields[9:] = samples
                        if reason == "Swapped genotypes.":
                            info += ";{0}{1}{2}{3}{4}{5}".format("AA=", allele_from_bed, ";CS=", cs_from_bed, ";SwappedGenotype=", "True")
                        else:
                            info += ";{0}{1}{2}{3}{4}{5}".format("AA=", allele_from_bed, ";CS=", cs_from_bed, ";SwappedGenotype=", "False")
                        fields[7] = info
                        out.write('\t'.join(fields) + '\n')
                    else:
                        reason = "No AA."
                else:
                    reason = "No CS."
            with open(fd_file, 'a') as fd:
                fd.write("Reason: " + reason + "\n")

            line_number = line_number + 1


AA_data = read_bed_file(args.bed)
print("BED file loaded. Filtering sites...")
update_vcf_with_bed_info(args.vcf, *AA_data, output_file = args.out, fd_file = args.fd)

print("Finished.")
