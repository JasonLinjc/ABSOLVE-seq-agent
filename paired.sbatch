#!/bin/bash
#SBATCH --job-name=absolve_pair
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=04:00:00
set -eo pipefail
# Slurm executes a spool copy. A direct invocation can resolve the checkout from this file.
agent_dir=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
# Manuscript mode orchestrates two software environments from the Slurm host.
panel=0
read_mode=
previous=
for arg in "$@"; do
    case "$arg" in
        --manuscript-panel) panel=1 ;;
        --read-mode=*) read_mode=${arg#*=} ;;
    esac
    if [[ "$previous" == --read-mode ]]; then read_mode=$arg; fi
    previous=$arg
done
if [[ "$read_mode" == manuscript || ( "$panel" == 1 && -z "$read_mode" ) ]]; then
    exec "${ABSOLVE_HOST_PYTHON:-python3}" -u "$agent_dir/paired.py" "$@"
fi
if [[ -n "${ABSOLVE_IMAGE:-}" ]]; then
    export APPTAINERENV_ABSOLVE_DATA="${ABSOLVE_DATA:-$agent_dir/data}"
    export APPTAINERENV_ABSOLVE_RUNS="${ABSOLVE_RUNS:-$agent_dir/runs}"
    export APPTAINERENV_SLURM_JOB_ID="$SLURM_JOB_ID"
    export APPTAINERENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}"
    for variable in ABSOLVE_MANIFEST ABSOLVE_CACHE ABSOLVE_MAXI_MERGED ABSOLVE_MAXI_POOLS; do
        if [[ -n "${!variable:-}" ]]; then
            export "APPTAINERENV_${variable}=${!variable}"
        fi
    done
    exec apptainer exec --cleanenv --bind "$agent_dir" "$ABSOLVE_IMAGE" \
        python -u "$agent_dir/paired.py" "$@"
fi
exec python -u "$agent_dir/paired.py" "$@"
