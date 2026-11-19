#!/bin/bash
#SBATCH --job-name=absolve_pair
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=04:00:00
set -eo pipefail
if [[ -n "${ABSOLVE_IMAGE:-}" ]]; then
    export APPTAINERENV_ABSOLVE_DATA="${ABSOLVE_DATA:-$SLURM_SUBMIT_DIR/data}"
    export APPTAINERENV_ABSOLVE_RUNS="${ABSOLVE_RUNS:-$SLURM_SUBMIT_DIR/runs}"
    export APPTAINERENV_SLURM_JOB_ID="$SLURM_JOB_ID"
    export APPTAINERENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}"
    for variable in ABSOLVE_MANIFEST ABSOLVE_CACHE ABSOLVE_MAXI_MERGED ABSOLVE_MAXI_POOLS; do
        if [[ -n "${!variable:-}" ]]; then
            export "APPTAINERENV_${variable}=${!variable}"
        fi
    done
    exec apptainer exec --cleanenv --bind "$SLURM_SUBMIT_DIR" "$ABSOLVE_IMAGE" \
        python -u /opt/absolve/paired.py "$@"
fi
exec python -u "${SLURM_SUBMIT_DIR}/paired.py" "$@"
