#!/bin/bash
#SBATCH --job-name=absolve-pbc
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -eu
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
exec python3 -u validate_pbc_panel.py "$@"
