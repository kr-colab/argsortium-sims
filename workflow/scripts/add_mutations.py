"""
Add mutations to simulated tree sequences using msprime.sim_mutations.
This script is called by Snakemake with parameters injected via the snakemake object.
"""
import stdpopsim
import msprime
import tszip
import time


def main():
    # Get parameters from Snakemake
    input_trees = snakemake.input.trees
    contig_name = snakemake.wildcards.contig
    species_name = snakemake.params.species
    genetic_map = snakemake.params.genetic_map
    mutation_model = snakemake.params.mutation_model
    discrete_genome = snakemake.params.discrete_genome
    base_seed = snakemake.params.seed

    # Output files
    output_trees = snakemake.output.trees
    output_log = snakemake.output.log

    # Calculate chromosome-specific seed (offset by 1000 to differentiate from ancestry seed)
    contig_list = snakemake.config["contigs"]
    contig_index = contig_list.index(contig_name)
    seed = base_seed + contig_index + 1000

    # Start timing
    start_time = time.time()

    # Load the tree sequence from simulate rule (ancestry only, no mutations)
    ts = tszip.load(input_trees)

    # Get mutation rate from stdpopsim species
    species = stdpopsim.get_species(species_name)
    contig = species.get_contig(contig_name, genetic_map=genetic_map)
    mutation_rate = contig.mutation_rate

    # Parse mutation model
    # Supported models: HKY, JC69, GTR, BinaryMutationModel
    if mutation_model == "HKY":
        # HKY with default kappa (transition/transversion ratio) of 2.0
        model = msprime.HKY(kappa=2.0)
    elif mutation_model == "JC69":
        model = msprime.JC69()
    elif mutation_model == "GTR":
        # GTR with equal rates (essentially JC69)
        model = msprime.GTR(relative_rates=[1, 1, 1, 1, 1, 1])
    elif mutation_model == "BinaryMutationModel":
        model = msprime.BinaryMutationModel()
    else:
        raise ValueError(
            f"Unknown mutation model: {mutation_model}. "
            f"Supported models: HKY, JC69, GTR, BinaryMutationModel"
        )

    # Add mutations to the tree sequence
    ts_mutated = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        model=model,
        discrete_genome=discrete_genome,
        random_seed=seed,
    )

    # Save mutated tree sequence
    tszip.compress(ts_mutated, output_trees)

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Write log file
    with open(output_log, "w") as f:
        f.write(f"Mutation Log\n")
        f.write(f"=" * 50 + "\n")
        f.write(f"Species: {species_name}\n")
        f.write(f"Contig: {contig_name}\n")
        f.write(f"Contig length: {contig.length:,} bp\n")
        f.write(f"Mutation rate: {mutation_rate:.2e} per bp per generation\n")
        f.write(f"Mutation model: {mutation_model}\n")
        f.write(f"Discrete genome: {discrete_genome}\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"\n")
        f.write(f"Input tree sequence:\n")
        f.write(f"  Number of trees: {ts.num_trees:,}\n")
        f.write(f"  Number of nodes: {ts.num_nodes:,}\n")
        f.write(f"  Number of samples: {ts.num_samples}\n")
        f.write(f"  Number of mutations: {ts.num_mutations:,}\n")
        f.write(f"\n")
        f.write(f"Output tree sequence (with mutations):\n")
        f.write(f"  Number of trees: {ts_mutated.num_trees:,}\n")
        f.write(f"  Number of nodes: {ts_mutated.num_nodes:,}\n")
        f.write(f"  Number of samples: {ts_mutated.num_samples}\n")
        f.write(f"  Number of mutations: {ts_mutated.num_mutations:,}\n")
        f.write(f"  Number of sites: {ts_mutated.num_sites:,}\n")
        f.write(f"\n")
        f.write(f"Timing:\n")
        f.write(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        f.write(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        f.write(f"  Elapsed time: {elapsed_time:.2f} seconds\n")


if __name__ == "__main__":
    main()
