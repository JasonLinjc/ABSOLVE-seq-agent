# Run D3734/T1 from the original FASTQs

This walkthrough runs the original D3734/T1 NoEP (LV_S3) and HiFi (LV_S1) libraries through all 79 targets. Each output row includes both indel rates and their paired difference (HiFi minus NoEP). The manuscript pBC rule and cross-target filter are the defaults.

## 1. Get the agent and create an environment

On a Linux Slurm cluster, choose a shared working directory visible to the login and compute nodes. The login node needs access to GitHub and the Conda channels. Use a new checkout and an empty environment prefix.

```bash
git clone https://github.com/JasonLinjc/ABSOLVE-seq-agent.git
cd ABSOLVE-seq-agent
unset ABSOLVE_DATA ABSOLVE_MANIFEST ABSOLVE_CACHE PYTHONPATH
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source "$(conda info --base)/etc/profile.d/conda.sh"
```

For a shared environment, run the following commands. For a node-local environment on Luminary, go directly to [Node-local installation](#node-local-installation).

```bash
env_prefix="$PWD/../absolve-pair-env"
mamba create --no-rc -y -p "$env_prefix" -f environment.yml
conda activate "$env_prefix"
python --version
CRISPResso --version
```

CRISPResso2 should report version 2.3.4. For the exact Linux x86-64 package builds used during acceptance, substitute `--file docs/environment-linux-64.lock` for `-f environment.yml`.

## 2. Build the published reference index

```bash
mkdir -p data
git init -q data/upstream
git -C data/upstream remote add origin https://github.com/pellinlab/ABSOLVE-seq.git
git -C data/upstream fetch --depth=1 origin 271355dfe40cf46b952e82b88dd4521f23f9feb4
git -C data/upstream checkout --detach FETCH_HEAD
python setup_references.py
```

The setup script checks the pinned sample and target tables and the SHA-256 of the published pBC annotation. It writes a fresh index in `data/references/`.

## 3. Identify the four input files

Set `FASTQ_ROOT` to the directory containing the two condition folders. Keep the original filenames. Each pattern below must match one file.

| Condition | Forward read | Reverse read |
| --- | --- | --- |
| NoEP | `NoEP/LV_S3_CKDL*_1.fq.gz` | `NoEP/LV_S3_CKDL*_2.fq.gz` |
| HiFi | `HiFi/LV_S1_CKDL*_1.fq.gz` | `HiFi/LV_S1_CKDL*_2.fq.gz` |

```bash
export FASTQ_ROOT=/path/to/original/fastqs
python examples/d3734_t1_manifest.py "$FASTQ_ROOT" data/sample_manifest.tsv
python paired.py --list
```

The example hashes the four compressed files and writes a two-row manifest. Compare its four SHA-256 values with your input inventory. Pair discovery should list `D3734`, `T1`, `NoEP_pool_D3734_T1` and `HiFi_pool_D3734_T1` together.

## 4. Submit the full pair

Submit from the agent repository root with the new environment active. Set the partition to one available on your cluster. The following command uses the Luminary CPU partitions.

```bash
run_dir="$PWD/runs/acceptance"
job_id=$(sbatch --parsable --partition=cpu1,cpu2 paired.sbatch \
    --donor D3734 --replicate T1 --fresh --output "$run_dir")
squeue -j "$job_id"
```

The job verifies both mates for each library, prepares all input pairs, and runs CRISPResso2 per library and target. Check for `COMPLETED` and exit code `0:0` before verifying results. Use a new output directory for each run.

```bash
sacct -j "$job_id" --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

## 5. Verify the result

Check both the table layout and its content. The SHA-256 below is the fingerprint of the independently verified D3734/T1 table for the original four FASTQs, the manuscript defaults and the v0.1.1 export format.

```bash
table="$run_dir/results/paired_target_rates.tsv"
python - "$table" <<'PY'
import csv
import sys
from collections import Counter

with open(sys.argv[1]) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

strategies = ("Unfiltered", "Dedup", "Dedud", "Dedud+Dedup")
keys = {(r["target"], r["strategy"]) for r in rows}
assert len(rows) == len(keys) == 316
assert len({t for t, _ in keys}) == 79
assert Counter(s for _, s in keys) == {s: 79 for s in strategies}
assert all((r["donor"], r["replicate"], r["maxi_method"], r["target_filter"])
           == ("D3734", "T1", "manuscript", "manuscript") for r in rows)
print("79 targets and 316 unique target/strategy rows.")
PY
printf '%s  %s\n' 0b0465ecc5aa244cc2a10654e153dbe8467322cbedd4f5bfc72806058018980b "$table" | sha256sum -c -
```

Read `NoEP_rate_pct` and `HiFi_rate_pct` as percentages. The `delta_HiFi_minus_NoEP_pp` column subtracts NoEP from HiFi and is measured in percentage points. A zero denominator leaves the corresponding rate and difference blank.

### Verified on Luminary

The fresh run used analysis release `v0.1.1` (`526d368`), CRISPResso2 2.3.4, six allocated CPUs and 18 GiB of requested memory. The paired analysis took 9 min 3 s. The run directory occupied 705 MiB. All 316 rows matched the archived paired table byte for byte. An independent calculation from the exported counts and denominators also reproduced both rates and every paired difference to within 1e-12.

| Library | Original read pairs | Pairs after preparation | Fresh target alignments |
| --- | ---: | ---: | ---: |
| NoEP | 10,084,551 | 5,971,199 | 79 |
| HiFi | 9,868,448 | 6,979,407 | 79 |

## Node-local installation

On Luminary, linking the environment on the shared filesystem was slow. The acceptance run installed the Linux packages on a compute node's local `/tmp` and ran the analysis in that allocation. The shared working directory held the FASTQs, reference index and final output.

After cloning the agent in step 1, run the Git commands in step 2 on the login node to fetch the upstream checkout. Defer `python setup_references.py` until the environment is active. Stage the locked archives in a shared cache and request an allocation from the agent repository root.

```bash
export CONDA_PKGS_DIRS="$PWD/data/conda-pkgs"
mkdir -p "$CONDA_PKGS_DIRS"
mamba create --no-rc --download-only -y -p "$PWD/data/download-env" \
    --file docs/environment-linux-64.lock
salloc --partition=cpu1,cpu2 --cpus-per-task=6 --mem=18G --time=04:00:00
srun --pty bash
```

In the allocated compute-node shell, load the same Mamba and Conda commands as on the login node. Convert the locked package URLs to local file URLs and install the environment. Keep this shell open through verification.

```bash
cd "$SLURM_SUBMIT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
export CONDA_PKGS_DIRS="$PWD/data/conda-pkgs"
scratch=$(mktemp -d /tmp/absolve.XXXXXXXX)
python3 - "$CONDA_PKGS_DIRS" "$scratch/local.lock" <<'PY'
from pathlib import Path
import sys
from urllib.parse import urlsplit

cache, output = map(Path, sys.argv[1:])
lines = ["@EXPLICIT"]
for entry in Path("docs/environment-linux-64.lock").read_text().splitlines():
    if entry.startswith("https://"):
        url = urlsplit(entry)
        package = cache / Path(url.path).name
        if not package.is_file():
            raise FileNotFoundError(package)
        lines.append(package.resolve().as_uri() + "#" + url.fragment)
output.write_text("\n".join(lines) + "\n")
PY
export MAMBA_ROOT_PREFIX="$scratch/mamba-root" CONDA_PKGS_DIRS="$scratch/pkgs"
mamba create --no-rc --offline -y -p "$scratch/env" --file "$scratch/local.lock"
conda activate "$scratch/env"
CRISPResso --version
```

Run `python setup_references.py` and the manifest commands in step 3. Then run the pair directly in the current allocation and verify the table using step 5.

```bash
run_dir="$PWD/runs/acceptance"
bash paired.sbatch --donor D3734 --replicate T1 --fresh --output "$run_dir"
```
