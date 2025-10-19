"""
Snakemake rules for adding mutations to simulated tree sequences.
"""

rule mutate:
    """
    Add mutations to simulated tree sequences using msprime.sim_mutations.
    Mutation rate is pulled from the stdpopsim Species object for each contig.
    """
    input:
        trees="{outdir}/{species}_{model}/{contig}/sim.trees"
    output:
        trees="{outdir}/{species}_{model}/{contig}/sim.mutated.trees",
        log="{outdir}/{species}_{model}/{contig}/sim.mutated.log"
    params:
        species=config["species_name"],
        genetic_map=config["genetic_map"],
        mutation_model=config.get("mutation_model", "HKY"),
        discrete_genome=config.get("discrete_genome", True),
        seed=config["base_seed"]
    log:
        "logs/{outdir}/{species}_{model}/{contig}.mutate.snakemake.log"
    script:
        "../scripts/add_mutations.py"
