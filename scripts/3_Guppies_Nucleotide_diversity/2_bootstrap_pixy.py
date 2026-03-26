import pandas as pd
import numpy as np
import subprocess
import os

# -------- USER SETTINGS --------
vcf_file = "Guppy_HQ.vcf.recode.vcf.gz"
pop_file = "pop.txt"
window_size = 1000000
n_boot = 100
n_cores = 1
# --------------------------------

np.random.seed(42)

# Load population file
pop_df = pd.read_csv(pop_file, sep="\t", header=None)
pop_df.columns = ["sample", "pop"]

results = []

for i in range(n_boot):

    print(f"Bootstrap replicate {i+1}/{n_boot}")

    boot_samples = []

    # Resample individuals WITH replacement per population
    for pop in pop_df["pop"].unique():
        sub = pop_df[pop_df["pop"] == pop]
        sampled = sub.sample(n=len(sub), replace=True)
        boot_samples.append(sampled)

    boot_df = pd.concat(boot_samples)

    # Write temporary population file
    boot_pop_file = f"boot_pop_{i}.txt"
    boot_df.to_csv(boot_pop_file, sep="\t", header=False, index=False)

    # Subset VCF
    sample_list = ",".join(boot_df["sample"].tolist())
    boot_vcf = f"boot_{i}.vcf.gz"

    subprocess.run([
        "bcftools", "view",
        "-s", sample_list,
        "-Oz",
        "-o", boot_vcf,
        vcf_file
    ], check=True)

    subprocess.run(["bcftools", "index", boot_vcf], check=True)

    # Run pixy
    out_folder = f"pixy_boot_{i}"
    subprocess.run([
        "pixy",
        "--stats", "pi",
        "--vcf", boot_vcf,
        "--populations", boot_pop_file,
        "--window_size", str(window_size),
        "--n_cores", str(n_cores),
        "--output_folder", out_folder,
        "--output_prefix", "boot"
    ], check=True)

    # Read pixy output
    pixy_out = os.path.join(out_folder, "boot_pi.txt")
    df_pixy = pd.read_csv(pixy_out, sep="\t")

    for pop in df_pixy["pop"].unique():
        df_pop = df_pixy[df_pixy["pop"] == pop]
        pi_val = df_pop["count_diffs"].sum() / df_pop["count_comparisons"].sum()
        results.append([i, pop, pi_val])

    # Cleanup (optional)
    os.remove(boot_pop_file)
    os.remove(boot_vcf)
    os.remove(boot_vcf + ".csi")

# Save results
results_df = pd.DataFrame(results, columns=["replicate", "pop", "pi"])
results_df.to_csv("bootstrap_pi_results.txt", sep="\t", index=False)

print("Bootstrap finished.")