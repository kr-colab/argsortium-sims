"""
Add mutations to simulated tree sequences using msprime.sim_mutations.
This script is called by Snakemake with parameters injected via the snakemake object.
"""
import tszip
import numpy as np
import gzip
import time


def main():
    # Get parameters from Snakemake
    input_trees = snakemake.input.trees
    contig_name = snakemake.wildcards.contig
    mispolarise = snakemake.params.mispolarise
    mask_missing = snakemake.params.mask_missing
    phasing_error = snakemake.params.phasing_error

    # Output files
    output_vcf = snakemake.output.vcf
    output_log = snakemake.output.log

    # Start timing
    start_time = time.time()

    # Load the tree sequence from simulate rule (ancestry only, no mutations)
    ts = tszip.load(input_trees)

    # Name individuals
    population_names = {p.id: p.metadata["name"] for p in ts.populations()}
    indiv_names = [f"{population_names[i.population]}_{i.id}" for i in ts.individuals()]

    # Add mispolarisation
    if mispolarise:
        assert False, "Not implemented"

    # Add masking
    if mask_missing:
        assert False, "Not implemented"

    # Add phasing error
    if phasing_error:
        assert False, "Not implemented"

    # Save VCF, positions are incremented by one relative to true tree sequence
    ts.write_vcf(
        gzip.open(output_vcf, "tw"),
        contig_id=contig_name,
        individual_names=indiv_names,
        position_transform=lambda x: np.array(x).astype(np.int64) + 1,
    )

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Write log file
    with open(output_log, "w") as f:
        f.write(f"VCF Log\n")
        f.write(f"=" * 50 + "\n")
        f.write(f"Input tree sequence:\n")
        f.write(f"  Number of trees: {ts.num_trees:,}\n")
        f.write(f"  Number of nodes: {ts.num_nodes:,}\n")
        f.write(f"  Number of samples: {ts.num_samples}\n")
        f.write(f"  Number of mutations: {ts.num_mutations:,}\n")
        f.write(f"\n")
        f.write(f"Timing:\n")
        f.write(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        f.write(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        f.write(f"  Elapsed time: {elapsed_time:.2f} seconds\n")


if __name__ == "__main__":
    main()
