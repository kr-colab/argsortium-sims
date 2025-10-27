"""
Snakemake rules for running population genetic simulations.
"""

rule simulate:
    """
    Run a single chromosome simulation using stdpopsim/msprime.
    Each chromosome runs independently and can be parallelized.
    """
    output:
        trees="{outdir}/{species}_{model}/{contig}/sim.trees",
        log="{outdir}/{species}_{model}/{contig}/sim.log",
        genmap="{outdir}/{species}_{model}/{contig}/sim.hapmap",
    params:
        species=config["species_name"],
        model=config["model_name"],
        samples=config["samples"],
        genetic_map=config["genetic_map"],
        add_outgroup=config["add_outgroup"],
        outgroup_time=config.get("outgroup_time", 215000),
        seed=config["base_seed"]
    log:
        "logs/{outdir}/{species}_{model}/{contig}.snakemake.log"
    script:
        "../scripts/run_simulation.py"
