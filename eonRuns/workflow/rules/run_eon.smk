# -*- mode:snakemake; -*-


rule do_minimization:
    input:
        config=lambda wildcards: f"resources/{config['compute_target']}/config_minim.ini",
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_pre_aligned.con",
        model=expand(
            f"{config['paths']['models']}/pet-mad-{{version}}.pt",
            version=config["pet_mad"]["version"],
        ),
    output:
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_minimized.con",
    shadow:
        "minimal"
    shell:
        """
        cp {input.config} config.ini
        cp {input.endpoint} pos.con
        cp -f {input.model} .
        eonclient
        cp min.con {output.endpoint}
        """


rule do_neb:
    input:
        config=lambda wildcards: f"resources/{config['compute_target']}/config_neb.ini",
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product.con",
        model=expand(
            f"{config['paths']['models']}/pet-mad-{{version}}.pt",
            version=config["pet_mad"]["version"],
        ),
    output:
        results_dat=f"{config['paths']['neb']}/{{system}}/results.dat",
        neb_con=f"{config['paths']['neb']}/{{system}}/neb.con",
        neb_dat=f"{config['paths']['neb']}/{{system}}/neb.dat",
    params:
        opath=f"{config['paths']['neb']}/{{system}}/",
    shell:
        """
        mkdir -p {params.opath}
        cp -f {input.model} {params.opath}/
        cp {input.config} {params.opath}/config.ini
        cp {input.reactant} {params.opath}/reactant.con
        cp {input.product} {params.opath}/product.con
        cd {params.opath}
        eonclient 2>&1 || true
        """
