# -*- mode:snakemake; -*-


rule plot_neb_path:
    input:
        con=f"{config['paths']['neb']}/neb.con",
        dat_files=expand(
            f"{config['paths']['neb']}/neb_{{num:03d}}.dat",
            num=range(config["eon"]["num_steps"]),
        ),
    output:
        plot1d_path=f"{config['paths']['plots']}/roneb_1D_path.png",
        plot1d_index=f"{config['paths']['plots']}/roneb_1D_index.png",
        plot1d_rmsd=f"{config['paths']['plots']}/roneb_1D_rmsd.png",
        plot2d_rmsd=f"{config['paths']['plots']}/roneb_2D_rmsd.png",
    params:
        plot_structures=config["plotting"]["plot_structures"],
        facecolor=config["plotting"]["facecolor"],
        figsize=config["plotting"]["figure"]["figsize"],
        dpi=config["plotting"]["figure"]["dpi"],
        fontsize_base=config["plotting"]["fonts"]["base"],
        title=config["plotting"]["title"],
        ipath=config["paths"]["neb"],
        zoom_ratio=config["plotting"]["figure"]["zoom_ratio"],
        dp_x=config["plotting"]["draw_product"]["x"],
        dp_y=config["plotting"]["draw_product"]["y"],
        dp_rad=config["plotting"]["draw_product"]["rad"],
        dr_x=config["plotting"]["draw_reactant"]["x"],
        dr_y=config["plotting"]["draw_reactant"]["y"],
        dr_rad=config["plotting"]["draw_reactant"]["rad"],
        ds_x=config["plotting"]["draw_saddle"]["x"],
        ds_y=config["plotting"]["draw_saddle"]["y"],
        ds_rad=config["plotting"]["draw_saddle"]["rad"],
        _aserot=config["plotting"]["aserot"],
        _cache=config["paths"]["cache"],
    shell:
        """
        mkdir -p {params._cache} &&
        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot1d_path} \
            --plot-type "profile" \
            --rc-mode "path" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} {params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} {params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --ase-rotation {params._aserot} &&


        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot1d_index} \
            --plot-type "profile" \
            --rc-mode "index" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} {params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} {params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --ase-rotation {params._aserot} &&


        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot1d_rmsd} \
            --plot-type "profile" \
            --rc-mode "rmsd" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} {params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} {params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --cache-file {params._cache}/1dcache.parquet \
            --ase-rotation {params._aserot} &&


        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot2d_rmsd} \
            --plot-type "landscape" \
            --rc-mode "rmsd" \
            --show-pts \
            --landscape-path "all" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --input-path-pattern "{params.ipath}/neb_path_*.con" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} -{params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} -{params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --cache-file {params._cache}/2dcache.parquet \
            --ase-rotation {params._aserot} \
            --rbf-smoothing 1e-4
        """
