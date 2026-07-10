# -*- mode:snakemake; -*-


def blog_system():
    return config["plotting"]["blog_system"]


def figure_dir():
    return f"{config['paths']['figures']}/blog"


def prefix():
    return config["plotting"]["blog_prefix"]


def system_from_landscape(wildcards):
    return {
        "cyclo": "system_100",
        "grignard": "11_grignard",
        "bicyclo": "bicyclobutane_05",
    }[wildcards.landscape]


def landscape_surface(wildcards):
    if system_from_landscape(wildcards) == "system_100":
        return config["plotting"]["surface_type_small"]
    return config["plotting"]["surface_type_large"]


def saddle_label(wildcards):
    # Short strip/legend labels (full citation lives in the post captions).
    return {
        "system_100": "ORCA",
        "11_grignard": "Birkholz",
        "bicyclobutane_05": "mlFSM",
    }[system_from_landscape(wildcards)]


def landscape_title(wildcards):
    return {
        "cyclo": "Cycloaddition",
        "grignard": "Grignard addition",
        "bicyclo": "Bicyclobutane opening",
    }[wildcards.landscape]


rule stack_info:
    output:
        f"{figure_dir()}/stack.json",
    shell:
        "python scripts/write_stack_info.py --output {output}"


rule plot_pipeline:
    output:
        f"{figure_dir()}/{prefix()}-pipeline.png",
    shell:
        "python scripts/plot_neb_diagnostics.py pipeline --output {output}"


rule plot_profile:
    input:
        con=f"{config['paths']['neb']}/{blog_system()}/neb.con",
        results=f"{config['paths']['neb']}/{blog_system()}/results.dat",
    output:
        f"{figure_dir()}/{prefix()}-1d-{{mode}}.png",
    wildcard_constraints:
        mode="path|index|rmsd",
    params:
        dat_pattern=f"{config['paths']['neb']}/{blog_system()}/neb_*.dat",
        cache=f"{config['paths']['cache']}/{blog_system()}/profile-{{mode}}.parquet",
        figsize=config["plotting"]["figsize_profile"],
        dpi=config["plotting"]["dpi"],
        fontsize_base=config["plotting"]["fontsize_base"],
        zoom_ratio=config["plotting"]["zoom_ratio"],
        facecolor=config["plotting"]["facecolor"],
        rotation=config["plotting"]["rotation"],
        plot_structures=config["plotting"]["plot_structures"],
        strip_renderer=config["plotting"].get("strip_renderer", "xyzrender"),
        perspective_tilt=config["plotting"].get("perspective_tilt", 8),
    shell:
        """
        mkdir -p $(dirname {params.cache}) &&
        python -m rgpycrumbs.cli eon plt-neb \
            --con-file {input.con} \
            --output-file {output} \
            --plot-type profile \
            --rc-mode {wildcards.mode} \
            --spline-method $([ "{wildcards.mode}" = "index" ] && echo spline || echo hermite) \
            --input-dat-pattern "{params.dat_pattern}" \
            --plot-structures {params.plot_structures} \
            --facecolor {params.facecolor} \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --cache-file {params.cache} \
            --rotation {params.rotation} \
            --strip-renderer {params.strip_renderer} \
            --perspective-tilt {params.perspective_tilt}
        """


rule plot_landscape:
    input:
        con=lambda wildcards: f"{config['paths']['neb']}/{system_from_landscape(wildcards)}/neb.con",
        results=lambda wildcards: f"{config['paths']['neb']}/{system_from_landscape(wildcards)}/results.dat",
        saddle=lambda wildcards: config["systems"][system_from_landscape(wildcards)]["saddle"],
    output:
        f"{figure_dir()}/{prefix()}-2d-{{landscape}}.png",
    wildcard_constraints:
        landscape="cyclo|grignard|bicyclo",
    params:
        dat_pattern=lambda wildcards: f"{config['paths']['neb']}/{system_from_landscape(wildcards)}/neb_*.dat",
        path_pattern=lambda wildcards: f"{config['paths']['neb']}/{system_from_landscape(wildcards)}/neb_path*.con",
        cache=lambda wildcards: f"{config['paths']['cache']}/{system_from_landscape(wildcards)}/landscape.parquet",
        figsize=config["plotting"]["figsize_landscape"],
        dpi=config["plotting"]["dpi"],
        fontsize_base=config["plotting"]["fontsize_base"],
        zoom_ratio=config["plotting"]["zoom_ratio"],
        facecolor=config["plotting"]["facecolor"],
        rotation=config["plotting"]["rotation"],
        plot_structures=config["plotting"]["plot_structures"],
        ira_kmax=config["plotting"]["ira_kmax"],
        surface=landscape_surface,
        label=saddle_label,
        title=landscape_title,
        strip_renderer=config["plotting"].get("strip_renderer", "xyzrender"),
        perspective_tilt=config["plotting"].get("perspective_tilt", 8),
        cmap_landscape=config["plotting"].get("cmap_landscape", "viridis"),
    shell:
        """
        mkdir -p $(dirname {params.cache}) &&
        python -m rgpycrumbs.cli eon plt-neb \
            --con-file {input.con} \
            --output-file {output} \
            --plot-type landscape \
            --project-path \
            --show-pts \
            --landscape-path all \
            --plot-structures {params.plot_structures} \
            --surface-type {params.surface} \
            --input-dat-pattern "{params.dat_pattern}" \
            --input-path-pattern "{params.path_pattern}" \
            --ira-kmax {params.ira_kmax} \
            --show-legend \
            --title "{params.title}" \
            --cmap-landscape {params.cmap_landscape} \
            --facecolor {params.facecolor} \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --cache-file {params.cache} \
            --rotation {params.rotation} \
            --strip-renderer {params.strip_renderer} \
            --perspective-tilt {params.perspective_tilt} \
            --additional-con {input.saddle} "{params.label}"
        """


rule readcon_summary:
    input:
        expand(
            f"{config['paths']['neb']}/{{system}}/neb.con",
            system=list(config["systems"].keys()),
        ),
    output:
        f"{figure_dir()}/readcon-v2-summary.csv",
    shell:
        "python scripts/plot_neb_diagnostics.py readcon-summary --output {output} {input}"


rule plot_convergence:
    input:
        expand(
            f"{config['paths']['neb']}/{{system}}/results.dat",
            system=list(config["systems"].keys()),
        ),
        expand(
            f"{config['paths']['neb']}/{{system}}/neb.con",
            system=list(config["systems"].keys()),
        ),
    output:
        f"{figure_dir()}/{prefix()}-convergence.png",
        f"{figure_dir()}/{prefix()}-force-heatmap.png",
        f"{figure_dir()}/{prefix()}-sampling-density.png",
    params:
        neb_root=config["paths"]["neb"],
        blog_system=blog_system(),
    shell:
        "python scripts/plot_neb_diagnostics.py diagnostics "
        "--neb-root {params.neb_root} "
        "--blog-system {params.blog_system} "
        "--convergence {output[0]} "
        "--force-heatmap {output[1]} "
        "--sampling-density {output[2]}"
