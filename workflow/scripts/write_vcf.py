"""
Add mutations to simulated tree sequences using msprime.sim_mutations.
This script is called by Snakemake with parameters injected via the snakemake object.
"""
import tszip
import numpy as np
import gzip
import textwrap
import time


def main():
    # Get parameters from Snakemake
    input_trees = snakemake.input.trees
    contig_name = snakemake.wildcards.contig
    add_outgroup = snakemake.params.add_outgroup
    increment_positions = snakemake.params.increment_positions
    mispolarise = snakemake.params.mispolarise
    mask_missing = snakemake.params.mask_missing
    add_phasing_error = snakemake.params.add_phasing_error

    # Output files
    output_vcf = snakemake.output.vcf
    output_log = snakemake.output.log
    output_outgroup = snakemake.output.outgroup_fasta
    output_ancestral = snakemake.output.ancestral_fasta

    # Start timing
    start_time = time.time()

    # Load the tree sequence from simulate rule (ancestry only, no mutations)
    ts = tszip.load(input_trees)

    # Store true ancestral states at each base
    ancestral_states = np.full(int(ts.sequence_length), "N")
    for s in ts.sites():
        ancestral_states[int(s.position)] = s.ancestral_state

    # Store outgroup state at each base and strip outgroup individual from tree sequence
    outgroup_states = np.full(int(ts.sequence_length), "N")
    if add_outgroup:
        samples = np.array(list(ts.samples()))
        outgroup_population, = [i for i,p in enumerate(ts.populations()) if p.metadata["name"] == "outgroup"]
        outgroup_samples = samples[ts.nodes_population[samples] == outgroup_population]
        for v in ts.variants(samples=outgroup_samples):
            outgroup_allele, = v.counts().most_common(1)
            outgroup_states[int(v.site.position)] = outgroup_allele[0]
        ts = ts.simplify(samples=[i for i in samples if not i in outgroup_samples])
    num_outgroup_mismatches = np.sum(outgroup_states != ancestral_states)

    # Name individuals
    population_names = {p.id: p.metadata["name"] for p in ts.populations()}
    indiv_names = [f"{population_names[i.population]}_{i.id}" for i in ts.individuals()]

    # Mask variants from VCF
    if mask_missing:
        assert False, "Not implemented"

    # Mispolarise VCF according to outgroup
    if add_outgroup and mispolarise:
        assert False, "Not implemented"

    # Add phasing error to VCF
    if add_phasing_error:
        assert False, "Not implemented"

    # Save VCF, positions are incremented by one relative to true tree sequence
    ts.write_vcf(
        gzip.open(output_vcf, "tw"),
        contig_id=contig_name,
        individual_names=indiv_names,
        position_transform=lambda x: np.array(x).astype(np.int64) + 1 if increment_positions else None,
    )

    # Save ancestral and outgroup states as gzipped fasta
    with gzip.open(output_ancestral, "tw") as f:
        f.write(f">{contig_name}\n" + "".join(ancestral_states) + "\n")
    with gzip.open(output_outgroup, "tw") as f:
        f.write(f">{contig_name}\n" + "".join(outgroup_states) + "\n")

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
        f.write(f"Ancestral states:\n")
        f.write(f"  Mismatches between outgroup and ancestral state: {num_outgroup_mismatches}\n")
        f.write(f"Timing:\n")
        f.write(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        f.write(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        f.write(f"  Elapsed time: {elapsed_time:.2f} seconds\n")


if __name__ == "__main__":
    main()
