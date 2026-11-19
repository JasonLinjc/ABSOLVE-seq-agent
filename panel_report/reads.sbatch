#!/bin/bash
#SBATCH --job-name=absolve_reads
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=08:00:00
set -euo pipefail
exec "${ABSOLVE_HOST_PYTHON:-python3}" "$1" "$2" "$3" "$SLURM_ARRAY_TASK_ID"
