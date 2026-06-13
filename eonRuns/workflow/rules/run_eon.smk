# -*- mode:snakemake; -*-
from rgpycrumbs.eon.helpers import write_eon_config
from pathlib import Path


rule do_minimization:
    input:
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_pre_aligned.con",
        model=expand(
            f"{config['paths']['models']}/{{model}}.pt",
            model=config["upet"]["model"],
        ),
    output:
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_minimized.con",
    params:
        device=config["compute_target"],
    shadow:
        "minimal"
    run:
        min_settings = {
            "Main": {"job": "minimization", "random_seed": 706253457},
            "Potential": {"potential": "metatomic"},
            "Metatomic": {
                "model_path": Path(input.model[0]).absolute(),
                "device": params.device,
            },
            "Optimizer": {
                "max_iterations": 2000,
                "opt_method": "lbfgs",
                "max_move": 0.1,
                "converged_force": 0.0514221,
            },
        }
        import shutil, subprocess

        write_eon_config(Path("."), min_settings)
        shutil.copy(input.endpoint, "pos.con")
        subprocess.run(["eonclient"], check=True)
        shutil.copy("min.con", output.endpoint)


rule do_neb:
    input:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product.con",
        model=expand(
            f"{config['paths']['models']}/{{model}}.pt",
            model=config["upet"]["model"],
        ),
    output:
        results_dat=f"{config['paths']['neb']}/{{system}}/results.dat",
        neb_con=f"{config['paths']['neb']}/{{system}}/neb.con",
        neb_dat=f"{config['paths']['neb']}/{{system}}/neb.dat",
    params:
        device=config["compute_target"],
        opath=f"{config['paths']['neb']}/{{system}}",
    run:
        neb_settings = {
            "Main": {"job": "nudged_elastic_band", "random_seed": 706253457},
            "Potential": {"potential": "metatomic"},
            "Metatomic": {
                "model_path": Path(input.model[0]).absolute(),
                "device": params.device,
            },
            "Nudged Elastic Band": {
                "images": 18,
                "energy_weighted": "true",
                "ew_ksp_min": 0.972,
                "ew_ksp_max": 9.72,
                "ew_trigger": 0.5,
                "initializer": "sidpp",
                "sidpp_growth_alpha": 0.33,
                "minimize_endpoints": "false",
                "climbing_image_method": "true",
                "climbing_image_converged_only": "true",
                "ci_after": 0.5,
                "ci_after_rel": 0.8,
                "ci_mmf": "true",
                "ci_mmf_after": 0.1,
                "ci_mmf_after_rel": 0.8,
                "ci_mmf_penalty_strength": 1.5,
                "ci_mmf_penalty_base": 0.4,
                "ci_mmf_angle": 0.9,
                "ci_mmf_nsteps": 1000,
                "ci_mmf_ci_stability_count": 5,
            },
            "Dimer": {
                "improved": "true",
                "opt_method": "cg",
                "remove_rotation": "false",
                "converged_angle": 10.0,
            },
            "Optimizer": {
                "max_iterations": 1000,
                "opt_method": "lbfgs",
                "max_move": 0.1,
                # marks and gomez use 10^-3 Ha / Bohr
                # 0.0514221
                "converged_force": 0.0514221,
            },
            "Debug": {"write_movies": "true"},
        }

        out_path = Path(params.opath)
        out_path.mkdir(parents=True, exist_ok=True)

        write_eon_config(out_path, neb_settings)

        import shutil, subprocess

        shutil.copy2(os.path.abspath(input.reactant), out_path / "reactant.con")
        shutil.copy2(os.path.abspath(input.product), out_path / "product.con")

        subprocess.run(["eonclient"], cwd=out_path, check=True)
