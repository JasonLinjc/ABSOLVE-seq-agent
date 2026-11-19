#!/bin/bash
#SBATCH --job-name=absolve_agent_pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -eo pipefail
python pilot.py --manifest "${ABSOLVE_MANIFEST:-${ABSOLVE_DATA:-data}/sample_manifest.tsv}" --output "runs/$SLURM_JOB_ID" "$@"
