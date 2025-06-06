import os
import glob
import csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class VariantCallingAnalysis:
    def __init__(self, results_dir="Results"):
        self.results_dir = os.path.abspath(results_dir)
        os.makedirs(self.results_dir, exist_ok=True)

        self.coverage_summaries = glob.glob(os.path.join(self.results_dir, "*.coverage_summary.txt"))
        self.depth_files = glob.glob(os.path.join(self.results_dir, "*.depth.txt"))
        self.vcf_files = glob.glob(os.path.join(self.results_dir, "*.vcf"))

    def open_file(self, location):
        return pd.read_csv(os.path.expanduser(location), sep="\t", header=None)

    def convert_all_to_csv(self):
        print("[INFO] Converting depth and coverage files to CSV...")
        for file in self.depth_files + self.coverage_summaries:
            try:
                df = self.open_file(file)
                if not df.empty:
                    csv_path = file.replace(".txt", ".csv")
                    df.to_csv(csv_path, index=False)
            except Exception as e:
                print(f"[WARNING] Could not convert {file} to CSV: {e}")

        print("[INFO] Converting VCF files to CSV...")
        for vcf_file in self.vcf_files:
            try:
                csv_file = vcf_file.replace(".vcf", ".csv")
                with open(vcf_file, "r") as vcf, open(csv_file, "w", newline="") as csv_out:
                    writer = csv.writer(csv_out)
                    for line in vcf:
                        if line.startswith("##"):
                            continue
                        elif line.startswith("#"):
                            header = line.strip("#").strip().split("\t")
                            writer.writerow(header)
                        else:
                            row = line.strip().split("\t")
                            writer.writerow(row)
            except Exception as e:
                print(f"[WARNING] Could not convert {vcf_file} to CSV: {e}")

        print("[DONE] Conversion complete.")

    def plot_coverage_summary(self):
        summary_data = []
        for file in self.coverage_summaries:
            sample = os.path.basename(file).split(".")[0]
            try:
                df = pd.read_csv(file, sep="\t")
                if df.shape[1] >= 2:
                    avg_depth = pd.to_numeric(df.iloc[0, 0], errors="coerce")
                    coverage = pd.to_numeric(df.iloc[0, 1], errors="coerce")
                    if pd.notnull(avg_depth) and pd.notnull(coverage):
                        summary_data.append([sample, avg_depth, coverage])
            except Exception as e:
                print(f"[WARNING] Skipped coverage summary for {sample}: {e}")
        if summary_data:
            df_summary = pd.DataFrame(summary_data, columns=["Sample", "Average Depth", "Coverage Breadth (%)"])
            df_summary.set_index("Sample").plot(kind="bar", figsize=(10, 6), title="Coverage Metrics per Sample")
            plt.ylabel("Value")
            plt.tight_layout()
            plt.savefig(os.path.join(self.results_dir, "coverage_metrics_per_sample.png"))
            print("[DONE] Saved: coverage_metrics_per_sample.png")
            return df_summary

    def plot_depth_distribution(self):
        for file in self.depth_files:
            sample = os.path.basename(file).split(".")[0]
            try:
                df = pd.read_csv(file, sep="\t", header=None, names=["chrom", "pos", "depth"])
                if df.empty:
                    continue
                plt.figure(figsize=(10, 5))
                plt.hist(df["depth"], bins=100, color='skyblue', edgecolor='black')
                plt.title(f"Depth Distribution - {sample}")
                plt.xlabel("Coverage Depth")
                plt.ylabel("Number of Bases")
                plt.tight_layout()
                plt.savefig(os.path.join(self.results_dir, f"{sample}_depth_distribution.png"))
                print(f"[DONE] Saved: {sample}_depth_distribution.png")
            except Exception as e:
                print(f"[WARNING] Skipped depth distribution for {sample}: {e}")

    def plot_variant_types(self):
        variant_stats = []
        for file in self.vcf_files:
            sample = os.path.basename(file).split(".")[0]
            try:
                vcf = pd.read_csv(file, comment="#", sep="\t", header=None)
                if vcf.shape[1] < 5:
                    continue
                vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
                vcf["Variant_Type"] = vcf.apply(
                    lambda row: "SNP" if len(str(row["REF"])) == 1 and len(str(row["ALT"])) == 1 else "Indel", axis=1
                )
                counts = vcf["Variant_Type"].value_counts()
                variant_stats.append({
                    "Sample": sample,
                    "SNPs": counts.get("SNP", 0),
                    "Indels": counts.get("Indel", 0)
                })
            except Exception as e:
                print(f"[WARNING] Skipped variant type analysis for {sample}: {e}")
        if variant_stats:
            df_variants = pd.DataFrame(variant_stats).set_index("Sample")
            df_variants.plot(kind="bar", stacked=True, figsize=(10, 6), title="SNPs vs Indels per Sample")
            plt.ylabel("Count")
            plt.tight_layout()
            plt.savefig(os.path.join(self.results_dir, "variant_types_per_sample.png"))
            print("[DONE] Saved: variant_types_per_sample.png")

    def plot_variant_density(self):
        for file in self.vcf_files:
            sample = os.path.basename(file).split(".")[0]
            try:
                vcf = pd.read_csv(file, comment="#", sep="\t", header=None)
                if vcf.shape[1] < 2:
                    continue
                vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
                vcf["BIN"] = (vcf["POS"] // 10000) * 10000
                density = vcf.groupby(["CHROM", "BIN"]).size().reset_index(name="Variant_Count")
                plt.figure(figsize=(12, 4))
                sns.lineplot(data=density, x="BIN", y="Variant_Count", hue="CHROM")
                plt.title(f"Variant Density per 10kb Window - {sample}")
                plt.xlabel("Genomic Position (10kb bins)")
                plt.ylabel("Variant Count")
                plt.tight_layout()
                plt.savefig(os.path.join(self.results_dir, f"{sample}_variant_density.png"))
                print(f"[DONE] Saved: {sample}_variant_density.png")
            except Exception as e:
                print(f"[WARNING] Skipped variant density plot for {sample}: {e}")

    def run_all(self):
        self.convert_all_to_csv()
        self.plot_coverage_summary()
        self.plot_depth_distribution()
        self.plot_variant_types()
        self.plot_variant_density()
        print("[✓] All analysis complete.")


if __name__ == "__main__":
    import sys
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "Results"
    analysis = VariantCallingAnalysis(results_dir=results_dir)
    analysis.run_all()