# -*- mode:snakemake; -*-

import ase.io

rule prepare_endpoints_pre:
    input:
        reactant=lambda wildcards: config["systems"][wildcards.system]["reactant"],
        product=lambda wildcards: config["systems"][wildcards.system]["product"],
    output:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant_pre_aligned.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product_pre_aligned.con",
    run:
        try:
            reactant_atm = ase.io.read(input.reactant)
            product_atm = ase.io.read(input.product)
        except FileNotFoundError as e:
            raise Exception(f"Endpoint configurations not found: {e}")

        # Define unit cell and center coordinates
        for atm in [reactant_atm, product_atm]:
            atm.set_cell([25, 25, 25])
            atm.center()

        ase.io.write(output.reactant, reactant_atm)
        ase.io.write(output.product, product_atm)

rule prepare_endpoints_post:
    input:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant_minimized.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product_minimized.con",
    output:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product.con",
    run:
        try:
            reactant_atm = ase.io.read(input.reactant)
            product_atm = ase.io.read(input.product)
        except FileNotFoundError as e:
            raise Exception(f"Minimized endpoints not found: {e}")

        # Define unit cell and center coordinates
        for atm in [reactant_atm, product_atm]:
            atm.set_cell([25, 25, 25])
            atm.center()

        ase.io.write(output.reactant, reactant_atm)
        ase.io.write(output.product, product_atm)
