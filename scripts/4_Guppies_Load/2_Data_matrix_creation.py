import argparse

parser = argparse.ArgumentParser(description="Data MATRIX creation. This script takes all info on variants. "
                                             "Python v3.10")
parser.add_argument('-v', '--vcf', help='Input vcf file. ', type=str, required=True)
parser.add_argument('-o', '--out', help='Output prefix.', type=str, required=True)
parser.add_argument('-p', '--pop', help='Population name.', type=str, required=True)
args = parser.parse_args()

print("Hello!")
print("We are running 'Data MATRIX creation' script")
print("Input File: ", args.vcf)
print("Output Prefix: ", args.out)
print("Working...")

chromosomes = {
    "name": ["NC_024331.1", "NC_024332.1", "NC_024333.1", "NC_024334.1", "NC_024335.1", "NC_024336.1", "NC_024337.1",
             "NC_024338.1", "NC_024339.1", "NC_024340.1", "NC_024341.1", "NC_024342.1", "NC_024343.1", "NC_024344.1",
             "NC_024345.1", "NC_024346.1", "NC_024347.1", "NC_024348.1", "NC_024349.1", "NC_024350.1", "NC_024351.1",
             "NC_024352.1", "NC_024353.1"],
    "end_1": [6823135, 9257309, 7053088, 6299440, 6781749, 6305835, 6282673, 5589281, 6823559, 6563959, 5775112,
              5287915,
              6704839, 5667792, 6128943, 6639948, 6157602, 4405330, 5694147, 5277088, 5154768, 5049758, 3616919],
    "end_2": [27292542, 37029235, 28212354, 25197759, 27126995, 25223339, 25130691, 22357124, 27294238, 26255838,
              23100446, 21151659, 26819358, 22671168, 24515770, 26559794, 24630407, 17621321, 22776590, 21108354,
              20619073, 20199032, 14467677],
    "end_1_region": ["centromere", "centromere", "centromere", "telomere", "telomere", "centromere",
                     "centromere", "centromere", "centromere", "centromere", "centromere", "centromere",
                     "centromere", "telomere", "centromere", "telomere", "centromere", "centromere",
                     "telomere", "centromere", "telomere", "centromere", "telomere"],
    "end_2_region": ["telomere", "telomere", "telomere", "centromere", "centromere", "telomere",
                     "telomere", "telomere", "telomere", "telomere", "telomere", "telomere",
                     "telomere", "centromere", "telomere", "centromere", "telomere", "telomere",
                     "centromere", "telomere", "centromere", "telomere", "centromere"]
}


def check_variant_region(chrom, pos):
    if chrom in chromosomes["name"]:
        index = chromosomes["name"].index(chrom)
        start_position = chromosomes["end_1"][index]
        end_position = chromosomes["end_2"][index]
        start_region = chromosomes["end_1_region"][index]
        end_region = chromosomes["end_2_region"][index]

        if pos <= start_position:
            return start_region
        elif pos >= end_position:
            return end_region
        else:
            return "middle"

    else:
        return "NA"


with open(args.vcf, 'r') as vcf, open(args.out, 'w') as out_aa:
    for variant in vcf:
        if variant.startswith('#'):
            continue
        else:
            fields = variant.strip().split('\t')
            chrom = fields[0].strip()
            pos = int(fields[1].strip())
            region = check_variant_region(chrom, pos)
            info = fields[7].strip().split(";")
            for pair in info:
                key_value = pair.split('=')
                if key_value[0] == 'CS':
                    cs = float(key_value[1])
                    break
            site = fields[7].strip().split("|")[1]
            effect = fields[7].strip().split("|")[2]
            genotypes = fields[9:]

            sample_index = 1
            for sample in genotypes:
                genotype = sample.split(':')[0]
                if genotype == "0/0" or genotype == "0|0":
                    count = 0
                elif genotype == "0/1" or genotype == "0|1":
                    count = 1
                elif genotype == "1/1" or genotype == "1|1":
                    count = 2
                else:
                    count = 0

                sample_file = f'sample{sample_index}_{args.pop}.out'
                with open(sample_file, 'a') as out:
                    print("Sample", sample_index, sep="", file=out, end="\t")
                    print(chrom, sep="", file=out, end="\t")
                    print(region, sep="", file=out, end="\t")
                    print(cs, sep="", file=out, end="\t")
                    if cs > 2.:
                        print("deletrious", file=out, end="\t")
                    else:
                        print("non-deletrious", file=out, end="\t")
                    print(site, file=out, end="\t")
                    print(effect, file=out, end="\t")
                    print(count, file=out, end="\t")
                    print(args.pop, file=out, end="\n")
                sample_index = sample_index + 1

# The final matrix will be huge. It will have rows with each derived alleles, with columns: Sample Chromosome Region CS Site_Load Site_Effect Allele_Count

print("Finished :)")
