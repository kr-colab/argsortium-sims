"""
Snakemake rules for writing to compressed VCF
"""

rule vcf:
    """
    Write tree sequence to VCF, applying perturbations such as mispolarisation,
    masking of variants, or phasing error.
    """
    input:
        trees="{outdir}/{species}_{model}/{contig}/sim.mutated.trees"
    output:
        vcf="{outdir}/{species}_{model}/{contig}/sim.mutated.vcf.gz",
        ancestral_fasta="{outdir}/{species}_{model}/{contig}/sim.mutated.ancestral.fa.gz",
        outgroup_fasta="{outdir}/{species}_{model}/{contig}/sim.mutated.outgroup.fa.gz",
        log="{outdir}/{species}_{model}/{contig}/sim.mutated.vcf.log",
    log:
        "logs/{outdir}/{species}_{model}/{contig}.vcf.snakemake.log"
    params:
        add_outgroup = config["add_outgroup"],
        increment_positions = config.get("increment_vcf_positions", True),
        mispolarise = config.get("mispolarise", False),  # FIXME: placeholder
        mask_missing = config.get("mask_missing", False),  # FIXME: placeholder
        add_phasing_error = config.get("add_phasing_error", False),  # FIXME: placeholder
    script:
        "../scripts/write_vcf.py"
