"""
Snakemake rules for adding mutations to simulated tree sequences.
"""

rule mutate:
    """
    Add mutations to simulated tree sequences using msprime.sim_mutations.
    Mutation rate is pulled from the stdpopsim Species object for each contig.
    """
    input:
        trees="{outdir}/{species}_{model}/{contig}_{start}_{end}/{pop_str}/sim_seed{seed}.init.trees"
    output:
        trees="{outdir}/{species}_{model}/{contig}_{start}_{end}/{pop_str}/sim_seed{seed}.mutated.trees",
        log="{outdir}/{species}_{model}/{contig}_{start}_{end}/{pop_str}/sim_seed{seed}.mutated.log"
    params:
        species=config["species_name"],
        genetic_map=config["genetic_map"],
        mutation_model=config.get("mutation_model", "HKY"),
        discrete_genome=config.get("discrete_genome", True)
    log:
        "logs/{outdir}/{species}_{model}/{contig}_{start}_{end}/{pop_str}/sim_seed{seed}.mutated.snakemake.log"
    script:
        "../scripts/add_mutations.py"

