import argparse

parser = argparse.ArgumentParser(description="FILTERING VCF COVERAGE - Here we filter VCF files by different sample \
min and max DP. Samples which failed the filtration will be transformed into missing genotypes. Python  v3.10")
parser.add_argument('-v', '--vcf', help='Input vcf file. Must be uncompressed.', type=str, required=True)
parser.add_argument('-o', '--out', help='Output vcf file.', type=str, required=True)
parser.add_argument('-d', '--dp', help='Tresholds of DP. Must be given as file with two values per sample in one line\
, minimum first, maximum second.', type=str, required=True)
parser.add_argument('-g', '--gq', help='Treshold of GQ. Must be given as single number.', type=str, required=True)
parser.add_argument('-f', '--fd', help='Check that shows failed samples (those turned to missing).', type=str,
                    required=False)
args = parser.parse_args()

print("Hello!")
print("We are running FILTERING VCF COVERAGE script")
print("Input File: ", args.vcf)
print("Output File: ", args.out)
print("Working...")

minDP = []
maxDP = []

with open(args.dp, 'r') as file:
    DP = file.readlines()

GQ = int(args.gq)
delimiter = ":"

for i, indv in enumerate(DP, start=0):
    minDP.append(int(indv.split()[0]))
    maxDP.append(int(indv.split()[1]))

print(minDP)
print(maxDP)

with open(args.vcf, 'r') as file:
    for v, variant in enumerate(file, start=0):
        if variant.startswith('#'):
            with open(args.out, 'a') as f:
                print(variant, file=f, end="")
        else:
            spl = variant.split("\t")
            indv_count = 0
            if spl[6] == ".":
                with open(args.out, 'a') as f:
                    print("\t".join(spl[0:9]), file=f, end="\t")
            format_field = spl[8].split(":")
            dp_index = format_field.index("DP")
            if "RGQ" in format_field:
                gq_index = format_field.index("RGQ")
            else:
                gq_index = format_field.index("GQ")
            missing_check = True
            for s, sample in enumerate(spl[9:], start=0):
                individual = sample.split(":")
                if "AD" in format_field and len(individual) < 4:
                    missing_check = False
                    missing_reason = "truncated GQ"
                elif "AD" not in format_field and len(individual) < 3:
                    missing_check = False
                    missing_reason = "truncated GQ"
                else:
                    if individual[0] != "./.":
                        if individual[dp_index] != ".":
                            dp = int(individual[dp_index])
                        else:
                            dp = 0
                        if individual[gq_index] != ".":
                            gq = int(individual[gq_index])
                        else:
                            gq = 0
                    if gq >= GQ and maxDP[indv_count] >= dp >= minDP[indv_count]:
                        missing_check = True
                    else:
                        missing_check = False
                        missing_reason = "quality fail"
                if missing_check:
                    with open(args.out, 'a') as f:
                        if s == len(spl[9:]) - 1:
                            print(sample, file=f, end="")
                        else:
                            print(sample, file=f, end="\t")
                else:
                    with open(args.out, 'a') as f:
                        if s == len(spl[9:]) - 1:
                            print("./.:" + ":".join(individual[1:]), file=f, end="")
                        else:
                            print("./.:" + ":".join(individual[1:]), file=f, end="\t")
                        if args.fd is not None:
                            with open(args.fd, 'a') as x:
                                print("Individual nr ", indv_count+1, file=x, end="\t")
                                print(sample, file=x, end="\t")
                                print("Min/Max DP: ", minDP[indv_count], "/", maxDP[indv_count], file=x, end="\t")
                                print("Missing reason: ", missing_reason, file=x)
                indv_count = indv_count + 1

print("Finished.")

