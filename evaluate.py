import argparse
import os
import pysam
from python import mutation_signatures as ms

def main():
    parser = argparse.ArgumentParser(description="Process variant file and generate mutation profiles.")
    parser.add_argument("-v", "--var-path", required=True, help="Path to the variant VCF/SNV file")
    parser.add_argument("-r", "--ref-path", required=True, help="Path to the reference genome FASTA file")
    parser.add_argument("-s", "--sample-name", help="Sample name (default: inferred from variant file name)", default=None)
    parser.add_argument("-o", "--outdir", help="Output directory (default: same as variant file)", default=None)

    args = parser.parse_args()

    # Infer sample name if not provided
    sample_name = args.sample_name
    if sample_name is None:
        sample_name = os.path.basename(args.var_path)
        for ext in [".gz", ".vcf", ".bcf", ".snv"]:
            if sample_name.endswith(ext):
                sample_name = sample_name[: -len(ext)]

    # Infer output directory if not provided
    out_dir = args.outdir or os.path.dirname(args.var_path)

    ref_genome = pysam.FastaFile(args.ref_path)

    variants_96c = ms.variants_mut_profile(args.var_path, sample_name, ref_genome)
    mutaion_counts_96c = ms.create_96c_mutation_counts(variants_96c)

    variants_96c.write_csv(f"{out_dir}/{sample_name}_96c_mutation_profiles.tsv", separator="\t")
    print(f"Saved mutation profiles to: {out_dir}/{sample_name}_96c_mutation_profiles.tsv")
    
    mutaion_counts_96c.write_csv(f"{out_dir}/{sample_name}_96c_mutation_count.tsv", separator="\t")
    print(f"Saved 96 channnel mutation counts to: {out_dir}/{sample_name}_96c_mutation_count.tsv")
    
    ms.SBS96_plot(mutaion_counts_96c[:, 2].to_numpy(), label=sample_name, height=5, width=16, s=10, xticks_label=False, file=f"{out_dir}/{sample_name}_96c_mutation_signature.pdf")
    print(f"Saved 96 channnel mutation counts to: {out_dir}/{sample_name}_96c_mutation_signature.pdf")

if __name__ == "__main__":
    main()