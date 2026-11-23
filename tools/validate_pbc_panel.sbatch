#!/bin/bash
#SBATCH --job-name=absolve-pbc
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -eu
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
exec python3 -u -m tools.validate_pbc_panel "$@"
