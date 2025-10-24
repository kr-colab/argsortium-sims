"""
Population genetic simulation using stdpopsim and msprime.
This script is called by Snakemake with parameters injected via the snakemake object.
"""
import stdpopsim
import msprime
import time
import sys


def ratemap_to_hapmap(
    ratemap: msprime.RateMap, 
    contig_name: str, 
    missing_as_zero: bool = False,
) -> str:
    """
    Write a recombination rate map into hapmap format.
    """
    physical_position = ratemap.position.astype(np.int64)
    scaled_rate = ratemap.rate * 1e8
    map_position = ratemap.get_cumulative_mass(physical_position) * 100
    hapmap = ["Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)"]
    if missing_as_zero: 
        scaled_rate[np.isnan(scaled_rate)] = 0.0
    if np.isnan(scaled_rate[-1]):  # handle trailing NaN
        scaled_rate[-1] = 0.0
    else:
        scaled_rate = np.append(scaled_rate, 0.0)
    for rate, pos, map in zip(scaled_rate, physical_position, map_position):
        if not np.isnan(rate): 
            hapmap.append(f"{contig_name}\t{pos}\t{rate:.10f}\t{map:.10f}")
    hapmap = "\n".join(hapmap) + "\n"
    return hapmap


def check_ratemap(ratemap: msprime.RateMap, path: str):
    ratemap_ck = msprime.RateMap.read_hapmap(path, map_col=3)
    assert np.allclose(ratemap.position, ratemap_ck.position)
    assert np.allclose(ratemap.rate, ratemap_ck.rate)
    ratemap_ck = msprime.RateMap.read_hapmap(path, rate_col=2)
    assert np.allclose(ratemap.position, ratemap_ck.position)
    assert np.allclose(ratemap.rate, ratemap_ck.rate)


def main():
    # Get parameters from Snakemake
    contig_name = snakemake.wildcards.contig
    species_name = snakemake.params.species
    model_name = snakemake.params.model
    sample_dict = snakemake.params.samples
    genetic_map = snakemake.params.genetic_map
    add_outgroup = snakemake.params.add_outgroup
    outgroup_time = snakemake.params.outgroup_time
    base_seed = snakemake.params.seed

    # Output files
    output_trees = snakemake.output.trees
    output_log = snakemake.output.log
    output_genmap = snakemake.output.genmap

    # Calculate chromosome-specific seed
    # Extract contig index from config's contig list
    contig_list = snakemake.config["contigs"]
    contig_index = contig_list.index(contig_name)
    seed = base_seed + contig_index

    # Start timing
    start_time = time.time()

    # Get stdpopsim model from catalog
    species = stdpopsim.get_species(species_name)
    model = species.get_demographic_model(model_name)
    contig = species.get_contig(contig_name, genetic_map=genetic_map)
    samples = model.get_sample_sets(sample_dict)

    # Save genetic map
    with open(output_genmap, "w") as f:
        f.write(ratemap_to_hapmap(contig.recombination_map))
    check_ratemap(contig.recombination_map, output_genmap)

    if add_outgroup:
        # Add an ancient sample from population 0 at specified time
        # representing the outgroup (e.g., chimp)
        samples.append(
            msprime.SampleSet(1, population=0, time=outgroup_time, ploidy=1)
        )

    # Simulate tree sequence with msprime using that model
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=model.model,
        sequence_length=contig.length,
        recombination_rate=contig.recombination_map,
        random_seed=seed,
    )

    if add_outgroup:
        # Dump the tree sequence tables
        new_tables = ts.dump_tables()

        # Get the outgroup node id
        outgroup_indices = []
        for i, node in enumerate(new_tables.nodes):
            if node.time == outgroup_time:
                outgroup_indices.append(i)

        # Replace outgroup node times with 0.0
        for ind in outgroup_indices:
            new_tables.nodes[ind] = new_tables.nodes[ind].replace(time=0.0)

        # Re-sort the tables and re-create the tree sequence
        new_tables.sort()
        ts = new_tables.tree_sequence()

    # Save tree sequence
    ts.dump(output_trees)

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Write log file
    with open(output_log, "w") as f:
        f.write(f"Simulation Log\n")
        f.write(f"=" * 50 + "\n")
        f.write(f"Species: {species_name}\n")
        f.write(f"Model: {model_name}\n")
        f.write(f"Contig: {contig_name}\n")
        f.write(f"Contig length: {contig.length:,} bp\n")
        f.write(f"Genetic map: {genetic_map}\n")
        f.write(f"Samples: {sample_dict}\n")
        f.write(f"Add outgroup: {add_outgroup}\n")
        if add_outgroup:
            f.write(f"Outgroup time: {outgroup_time:,} generations\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"\n")
        f.write(f"Results:\n")
        f.write(f"  Number of trees: {ts.num_trees:,}\n")
        f.write(f"  Number of nodes: {ts.num_nodes:,}\n")
        f.write(f"  Number of samples: {ts.num_samples}\n")
        f.write(f"\n")
        f.write(f"Timing:\n")
        f.write(f"  Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        f.write(f"  End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        f.write(f"  Elapsed time: {elapsed_time:.2f} seconds\n")


if __name__ == "__main__":
    main()
