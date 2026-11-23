#!/bin/bash
#SBATCH --job-name=absolve_pair
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=04:00:00
set -eo pipefail
# Slurm executes a spool copy. A direct invocation can resolve the checkout from this file.
agent_dir=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
exec "${ABSOLVE_HOST_PYTHON:-python3}" -u "$agent_dir/paired.py" "$@"
