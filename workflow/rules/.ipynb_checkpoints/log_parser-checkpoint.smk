"""
Snakemake rules for parsing the log files and writing out a csv params for downstream inference
"""
rule log_parser:
    input:
        init_log="{outdir}/{species}_{model}/{pop_str}/{contig}_{start}_{end}_sim_seed{seed}.init.log",
        mutations_log="{outdir}/{species}_{model}/{pop_str}/{contig}_{start}_{end}_sim_seed{seed}.mutated.log",
        vcf_log="{outdir}/{species}_{model}/{pop_str}/{contig}_{start}_{end}_sim_seed{seed}.mutated.vcf.log"
    output:
        params_file="{outdir}/{species}_{model}/{pop_str}/{contig}_{start}_{end}_sim_seed{seed}.mutated.params.csv"
    shell:
        """
        python scripts/parse_logs_argparse.py \
            --init-log {input.init_log} \
            --mutations-log {input.mutations_log} \
            --vcf-log {input.vcf_log} \
            --output {output.params_file}
        """