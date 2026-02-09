# -*- mode:snakemake; -*-

rule convert_ckpt_to_pt:
    output:
        protected(f"{config['paths']['models']}/{{model}}.pt"),
    shell:
        "mtt export lab-cosmo/upet models/{wildcards.model}.ckpt && mv *.pt {output}"
