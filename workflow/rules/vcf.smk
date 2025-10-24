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
        log="{outdir}/{species}_{model}/{contig}/sim.mutated.vcf.log",
    log:
        "logs/{outdir}/{species}_{model}/{contig}.vcf.snakemake.log"
    params:
        mispolarise = False,  # FIXME: placeholder
        mask_missing = False,  # FIXME: placeholder
        phasing_error = False,  # FIXME: placeholder
    script:
        "../scripts/write_vcf.py"
