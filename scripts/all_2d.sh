#!/usr/bin/env bash

KERNELS=("grid" "rbf" "grad_matern" "grad_imq" "matern" "imq" "grad_rq" "grad_se")

for KERNEL in "${KERNELS[@]}"; do
    echo "Processing landscape with kernel: ${KERNEL}"

    python -m rgpycrumbs.cli eon plt_neb \
        --con-file "neb.con" \
        --output-file "landscape_${KERNEL}.png" \
        --plot-type "landscape" \
        --rc-mode "rmsd" \
        --show-pts \
        --landscape-path "all" \
        --plot-structures "crit_points" \
        --facecolor "white" \
        --figsize 5.37 5.37 \
        --dpi 300 \
        --zoom-ratio 0.25 \
        --fontsize-base 12 \
        --title "" \
        --ase-rotation 0x,0y,0z \
        --ira-kmax 14 \
        --surface-type "${KERNEL}" \
        --show-legend \
        --additional-con ../../resources/ewNEB/configurations/system100-sp.xyz "ORCA"
done
