#!/usr/bin/env bash
# Generates the 2D landscape figures for the paper.
# Requires rgpycrumbs >= 0.5.0 and chemparseplot >= 1.2.0 with the
# Nystrom-accelerated gradient-enhanced GP surface type (grad_imq_ny).
# Run from the nebviz_repro root inside `pixi s -e cpu`.
set -euo pipefail

RESULTS="eonRuns/results/02_neb"
IMGDIR="imgs"
mkdir -p "${IMGDIR}"

COMMON_OPTS=(
    --plot-type landscape
    --project-path
    --show-pts
    --landscape-path all
    --plot-structures crit_points
    --facecolor white
    --figsize 5.37 5.37
    --dpi 300
    --zoom-ratio 0.25
    --fontsize-base 12
    --title ""
    --ase-rotation 0x,0y,0z
    --show-legend
)

# 1. Cycloaddition (ethylene + N2O) -- full GP, small dataset
echo "=== Cycloaddition ==="
(cd "${RESULTS}/system_100" && \
    python -m rgpycrumbs.cli --dev eon plt-neb \
        --con-file neb.con \
        --output-file "${IMGDIR}/jax_grad_imq_cyclo5.png" \
        "${COMMON_OPTS[@]}" \
        --ira-kmax 14 \
        --surface-type grad_imq \
        --additional-con ../../resources/ewNEB/configurations/system100-sp.xyz "ORCA")

# 2. Grignard rearrangement -- Nystrom, large dataset
echo "=== Grignard ==="
(cd "${RESULTS}/11_grignard" && \
    python -m rgpycrumbs.cli --dev eon plt-neb \
        --con-file neb.con \
        --output-file "${IMGDIR}/jax_grad_imq_grignard.png" \
        "${COMMON_OPTS[@]}" \
        --ira-kmax 14 \
        --surface-type grad_imq_ny \
        --additional-con ts.xyz "Birkholz")

# 3. Bicyclobutane ring opening -- Nystrom, large dataset
echo "=== Bicyclobutane ==="
(cd "${RESULTS}/bicyclobutane_05" && \
    python -m rgpycrumbs.cli --dev eon plt-neb \
        --con-file neb.con \
        --output-file "${IMGDIR}/jax_grad_imq_bicyclobutane.png" \
        "${COMMON_OPTS[@]}" \
        --ira-kmax 14 \
        --surface-type grad_imq_ny \
        --additional-con ts.xyz "mlFSM")

echo "Done. Figures in ${IMGDIR}/."
