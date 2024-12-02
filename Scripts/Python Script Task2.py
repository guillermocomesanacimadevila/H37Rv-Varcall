# =========== Call libraries =========== #

import os
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt

# =========== Open files and change column names =========== #

def open_file(location):
    file = open(os.path.expanduser(location), 'r')
    return pd.read_csv(file, sep="\t", header=None)

depth_file_1 = open_file("~/Desktop/(MSc) BIOINFORMATICS/MACHINE LEARNING RESEARCH PROJECT (1B)/task2/coverage_1kb_means.bed")
depth_file_2 = open_file("~/Desktop/(MSc) BIOINFORMATICS/MACHINE LEARNING RESEARCH PROJECT (1B)/task2/coverage_1kb_means_2.bed")
depth_file_3 = open_file("~/Desktop/(MSc) BIOINFORMATICS/MACHINE LEARNING RESEARCH PROJECT (1B)/task2/coverage_1kb_means_3.bed")

column_names = ["reference_genome", "start", "end", "coverage_mean"]
depth_file_1.columns = column_names
depth_file_2.columns = column_names
depth_file_3.columns = column_names
print(depth_file_1)
# print(depth_file_1.head())

# Ensure dimensional shape equality
print(np.shape(depth_file_1))
print(np.shape(depth_file_2))
print(np.shape(depth_file_3))

# =========== Max and Min Coverage point filter =========== #

print(round(depth_file_1["coverage_mean"].min(), 3)) # 25.306
print(round(depth_file_1["coverage_mean"].min(), 3)) # 25.306
print(round(depth_file_1["coverage_mean"].min(), 3)) # 25.306

print(depth_file_1["coverage_mean"].max()) # 106.097 #@#
print(depth_file_2["coverage_mean"].max()) # 62.035
print(depth_file_3["coverage_mean"].max()) # 63.129

# Average depth per 1,000bp
print(depth_file_1.coverage_mean.mean()) # 50.10580652560685
print(depth_file_2.coverage_mean.mean()) # 49.996837682942626
print(depth_file_3.coverage_mean.mean()) # 49.99690890871785

max_coverage_mean = depth_file_1["coverage_mean"].max()
max_row = depth_file_1.loc[depth_file_1["coverage_mean"] == max_coverage_mean]
start_pos = max_row["start"].values
end_pos = max_row["end"].values

sort_1 = sorted(depth_file_1["coverage_mean"]) [::-1]
sort_2 = sorted(depth_file_2["coverage_mean"]) [::-1]
sort_3 = sorted(depth_file_3["coverage_mean"]) [::-1]

print(f"Sample 1 values: {[value for value in sort_1 if value > depth_file_1.coverage_mean.mean()]}")
print(f"Sample 2 values: {[value for value in sort_2 if value > depth_file_2.coverage_mean.mean()]}")
print(f"Sample 3 values: {[value for value in sort_3 if value > depth_file_3.coverage_mean.mean()]}")

# =========== Starting and Ending position for coverage outlier =========== #

# print(start_pos, end_pos)
print(type(depth_file_1)) # Printing type just out of interest
print(f"Maximum coverage_mean (file1): {max_coverage_mean}, Start: {start_pos}, End: {end_pos}")


# =========== Look at position 3406 on .bam file (A) =========== #

# Deletion (when counting SNP distance not count deletions)
# Compare all 3 seq with reference genome and invidivually


# =========== Analysing .vcf files =========== #
# T2Q2: Find ways iin which the reference genome might be wrong

# Contamination of reference genome wtihin that area specifically (the same everywhere except there) or it´s evolved

def vcf_to_csv(vcf_file_path, csv_file_path):
    with open(vcf_file_path, "r") as vcf_file:
        with open(csv_file_path, "w", newline="") as csv_file:
            csv_writer = csv.writer(csv_file)
            for line in vcf_file:
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    header = line.strip("#").strip().split("\t")
                    csv_writer.writerow(header)
                else:
                    data = line.strip().split("\t")
                    csv_writer.writerow(data)
            return data

# =========== Calling Variants - .csv conversion =========== #

vcf_file_path_1 = "~/Desktop/variants_1.vcf"
csv_file_path_1 = "~/Desktop/variants_1.csv"

vcf_file_path_2 = "~/Desktop/variants_2.vcf"
csv_file_path_2 = "~/Desktop/variants_2.csv"

vcf_file_path_3 = "~/Desktop/variants_3.vcf"
csv_file_path_3 = "~/Desktop/variants_3.csv"

# =========== vcf file 1 =========== #

vcf_to_csv(os.path.expanduser(vcf_file_path_1), os.path.expanduser(csv_file_path_1))
vcf_csv_path_1 = os.path.expanduser(csv_file_path_1)
vcf_data_1 = pd.read_csv(vcf_csv_path_1)
print(vcf_data_1.head())
print(type(vcf_data_1)) # Check it´s changed into a pd data frame
print(np.shape(vcf_data_1))

# Identifying positions where "ALT" allele == indel
positions_allele = vcf_data_1["ALT"]
indels = [allele for allele in positions_allele if len(allele) > 1]
positions = vcf_data_1.loc[vcf_data_1["ALT"].isin(indels)]
indel_pos_dict = {allele: positions.loc[positions["ALT"] == allele, "POS"].values.tolist() for allele in indels}
print(f"dictionary 1: {indel_pos_dict}")

min_distance = [min(value) for key, value in indel_pos_dict.items()]
indel_counter = 0
for indel in indels:
    indel_counter += 1

# =========== Map the depth to the corresponding indel ===========
# indel_coverage_dict = {allele: positions.loc[positions["coverage_mean"] == allele, "POS"].values.tolist() for allele in indels}
# print(f"Indel coverage mean: {indel_coverage_dict}")
# =========== Map the depth to the corresponding indel ===========

indel_coverage_dict = {allele: positions.loc[positions["ALT"] == allele, "DP"].values.tolist() for allele in indels} ###
print(f"Dictionary: {indel_coverage_dict}") ###

print(indel_counter) # number of indels (42)
print(min_distance[0], min_distance[-1]) # 1245282 & 1254357

# =========== vcf file 2 =========== #

vcf_to_csv(os.path.expanduser(vcf_file_path_2), os.path.expanduser(csv_file_path_2))
vcf_csv_path_2 = os.path.expanduser(csv_file_path_2)
vcf_data_2 = pd.read_csv(vcf_csv_path_2)
# print(vcf_data_2.head())
# print(type(vcf_data_2))

# =========== vcf file 3 =========== #

vcf_to_csv(os.path.expanduser(vcf_file_path_3), os.path.expanduser(csv_file_path_3))
vcf_csv_path_3 = os.path.expanduser(csv_file_path_3)
vcf_data_3 = pd.read_csv(vcf_csv_path_3)
# print(vcf_data_3.head())
# print(type(vcf_data_3))

# =========== Plotting =========== #

# Remember that the depth files have already been converted into pd dataframes
def plot_coverage(depth_file, plot_file=""):

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(depth_file["start"], depth_file["coverage_mean"], marker='o', linestyle="", color="b")
    ax.set_xlabel("Position (start of 1kb bin)")
    ax.set_ylabel("Mean Depth")
    ax.set_title("Mean Depth of Coverage in 1kb Bins")
    fig.tight_layout()

    if plot_file:
        fig.savefig(plot_file)
        print(f"Plot saved to '{plot_file}'.")
    else:
        plt.show()

    return fig, ax

fig, ax = plot_coverage(depth_file_1) # Visual outliers
fig, ax = plot_coverage(depth_file_2)
fig, ax = plot_coverage(depth_file_3)