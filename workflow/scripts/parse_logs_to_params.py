"""
Parse logs to params file for downstream analysis
This script uses argparse for command-line argument parsing.
"""
import pandas as pd
import re
import argparse


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Parse simulation logs and extract parameters to CSV file")
    parser.add_argument("--init-log", required=True, help="Path to initialization log file")
    parser.add_argument("--mutations-log", required=True, help="Path to mutations log file")
    parser.add_argument("--vcf-log", required=True, help="Path to VCF log file")
    parser.add_argument("--output", "-o", required=True, help="Path to output params CSV file")
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    init_log_path = args.init_log
    mutations_log_path = args.mutations_log
    vcf_log_path = args.vcf_log
    params_file = args.output
    
    # Read files
    with open(init_log_path, "r") as f:
        init_text = f.read()
    
    with open(mutations_log_path, "r") as f:
        mutated_text = f.read()
    
    with open(vcf_log_path, "r") as f:
        vcf_text = f.read()
    
    # ---- Extract fields ----
    # Model
    model_match = re.search(r"Model:\s*([A-Za-z0-9_]+)", init_text)
    model = model_match.group(1) if model_match else None
    
    # Contig
    contig_match = re.search(r"Contig:\s*(\w+)", init_text)
    contig = contig_match.group(1) if contig_match else None
    
    # Simulated sequence length
    sim_len_match = re.search(r"Simulated sequence length:\s*([0-9]+)\s*bp", init_text)
    sim_len = sim_len_match.group(1) if sim_len_match else None
    
    # Recombination / genetic map
    recomb_map_match = re.search(r"Genetic map written to:\s*(.*)", init_text)
    recomb_map = recomb_map_match.group(1) if recomb_map_match else None
    
    # Samples dictionary
    samples_match = re.search(r"Samples:\s*(\{.*?\})", init_text)
    samples = samples_match.group(1) if samples_match else None
    
    # Number of nodes (same in init + mutated, so take from init)
    leaves_match = re.search(r"Number of samples:\s*([0-9]+)", init_text)
    leaves = int(leaves_match.group(1)) if leaves_match else None
    
    # Random seed take from init
    seed_match = re.search(r"Random seed:\s*([0-9.eE+-]+)", init_text)
    seed = int(seed_match.group(1)) if seed_match else None
    
    # Mutation rate (μ) — only present in the mutated log
    mu_match = re.search(r"Mutation rate:\s*([0-9.eE+-]+)", mutated_text)
    mu = float(mu_match.group(1)) if mu_match else None
    
    # Ancestral fasta — only present in the vcf log
    ancestral_match = re.search(r"ancestral fasta written to:\s*(.*)", vcf_text)
    ancestral_fasta = ancestral_match.group(1).strip() if ancestral_match else None

    # VCF file path — only present in the vcf log
    vcf_match = re.search(r"VCF written:\s*(.*)", vcf_text)
    vcf_file = vcf_match.group(1).strip() if vcf_match else None
    
    # Build dataframe row
    row = {
        "vcf_file" : vcf_file,
        "mu": mu,
        "contig": contig,
        "samples": samples,
        "simulated_length": sim_len,
        "recomb_map": recomb_map,
        "leaves": leaves,
        "model": model,
        "seed": seed,
        "ancestral_fasta": ancestral_fasta,
    }
    
    pd.DataFrame([row]).to_csv(params_file, header=True, index=False)
    print(f"Parameters successfully written to {params_file}")


if __name__ == "__main__":
    main()
