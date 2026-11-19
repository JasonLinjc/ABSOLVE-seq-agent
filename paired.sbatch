#!/bin/bash
#SBATCH --job-name=absolve_pair
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=04:00:00
set -eo pipefail
exec python -u "${SLURM_SUBMIT_DIR}/paired.py" "$@"
