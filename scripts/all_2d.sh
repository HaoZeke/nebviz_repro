#!/usr/bin/env bash
# Generates the blog figures through the Snakemake workflow.
# Run from the repository root inside a pixi environment.
set -euo pipefail

cd eonRuns
snakemake --configfile config/general_config.yml -c4 all_figures
