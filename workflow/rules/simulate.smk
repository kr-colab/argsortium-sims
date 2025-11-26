rule simulate:
    """
    Run a single chromosome simulation using stdpopsim/msprime.
    """
    output:
        trees=(
            "{outdir}/{species}_{model}/{contig}_{start}_{end}/"
            "{pop_str}/sim_seed{seed}.init.trees"
        ),
        log=(
            "{outdir}/{species}_{model}/{contig}_{start}_{end}/"
            "{pop_str}/sim_seed{seed}.init.log"
        ),
        genmap=(
            "{outdir}/{species}_{model}/{contig}_{start}_{end}/"
            "{pop_str}/sim_seed{seed}.hapmap"
        )
    params:
        species=config["species_name"],
        model=config["model_name"],
        samples=config["samples"],
        genetic_map=config["genetic_map"],
        add_outgroup=config["add_outgroup"],
        outgroup_time=config.get("outgroup_time", 215000)
    script:
        "../scripts/run_simulation.py"
