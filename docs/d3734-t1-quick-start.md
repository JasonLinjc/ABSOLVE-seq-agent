# Run D3734/T1 from paired FASTQs

This walkthrough runs the original D3734/T1 NoEP (LV_S3) and HiFi (LV_S1) libraries through all 79 targets. Each output row includes both indel rates and their paired difference (HiFi minus NoEP). The manuscript pBC rule and cross-target filter are the defaults.

For a prebuilt environment, use the [container instructions](container.md). The published SIF reproduced the paired rate table on Luminary with the full public ENA downloads. Use the current checkout for the report, power estimates and figure.

## 1. Get the agent and create an environment

Use a local working directory on a workstation or a shared working directory on a Slurm cluster. Download the inputs on a machine with access to GitHub, ENA and the Conda channels. Use a new checkout and an empty environment prefix.

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

## 3. Get the four input files

### Public test data from SRA

ABSOLVE-seq study SRP722226 (BioProject PRJNA1503690) includes NoEP, HiFi Cas9 and wild-type Cas9 (WT) pooled reporter libraries from donor D3734, replicate T1. The paired quick start uses NoEP and HiFi to calculate HiFi minus NoEP across 79 targets. These four compressed ENA files total 1.58 GB. [NCBI study](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP722226) · [ENA file report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP722226&result=read_run&fields=run_accession,library_name,fastq_ftp,fastq_md5,fastq_bytes)

| Library | SRA run | Read pairs |
| --- | --- | ---: |
| `NoEP_pool_D3734_T1` | [SRR39887646](https://www.ebi.ac.uk/ena/browser/view/SRR39887646) | 10,084,551 |
| `HiFi_pool_D3734_T1` | [SRR39887695](https://www.ebi.ac.uk/ena/browser/view/SRR39887695) | 9,868,448 |
| `WT_pool_D3734_T1` | [SRR39887693](https://www.ebi.ac.uk/ena/browser/view/SRR39887693) | 19,418,487 |

The 26-library manuscript analysis includes WT and estimates WT minus NoEP across donors and replicates. See the [manuscript comparison](../README.md#manuscript-comparison) for the full input requirements.

For the paired quick start, download both mates of the NoEP and HiFi runs from the repository root on a networked node. Keep their ENA filenames. The manifest helper checks each compressed file against its published MD5 and records its SHA-256 for subsequent runs. [Download URLs and checksums](../examples/d3734_t1_ena.tsv)

```bash
export FASTQ_ROOT="$PWD/data/sra/SRP722226"
mkdir -p "$FASTQ_ROOT"
awk -F '\t' 'NR > 1 { print $4; print $6 }' examples/d3734_t1_ena.tsv |
    (cd "$FASTQ_ROOT" && xargs -n 1 -P 4 curl -fL --retry 3 --retry-all-errors -C - -O)
python examples/d3734_t1_manifest.py --ena "$FASTQ_ROOT" data/sample_manifest.tsv
python paired.py --list
```

If a connection drops partway through a file, keep the partial file. From `$FASTQ_ROOT`, rerun `curl -fL --retry 3 --retry-all-errors -C - -O "$mate_url"` with that mate’s URL from the download table, then run the manifest command. The helper verifies all four MD5 values before writing the manifest. Pair discovery should show `D3734 T1 NoEP_pool_D3734_T1 HiFi_pool_D3734_T1`. Continue with step 4.

### Existing laboratory FASTQs

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

The example hashes the four compressed files and writes a two-row manifest. Compare its four SHA-256 values with your laboratory inventory.

## 4. Run the full pair

The local command returns after it writes the report, tables, figure and Excel workbook. Use a new output directory.

```bash
run_dir="$PWD/runs/acceptance"
python paired.py --executor local --donor D3734 --replicate T1 --threads 4 --fresh --output "$run_dir"
python verify_results.py pair-fast "$run_dir/results"
python verify_results.py power "$run_dir/results"
```

### Slurm submission

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

## 5. Read the report and verify the result

Open `$run_dir/results/report.md` for QC, observed editing effects and estimated detection power. Its Dedud figure, PDF and Excel workbook are under `results/figures/`. The four-strategy power table has 632 rows (79 targets × 4 strategies × 2 effect sizes). [Power method](power.md)


Check both the table layout and its content. The original laboratory FASTQs and the four public ENA FASTQs produced the same 316-row table with the manuscript defaults. The checksum below applies to both input routes.

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

A clean-checkout run on 2026-11-21 tested the [one-message workflow](../README.md#start-with-one-message) at agent revision `681f213`. All four ENA FASTQs were downloaded on a networked workstation, transferred to Luminary, and verified from the cluster manifest builder. The reference index was built in the new checkout. Slurm job `21298427` used the published Apptainer image and `--fresh`, completed in 10 min 13 s (analysis only), and produced 79 new alignments per arm.

| Check | Result |
| --- | --- |
| Paired effects | All 316 rows matched the accepted SHA-256. Recalculation from the counts reproduced every arm rate and paired difference within 1e-12. |
| Power | All 632 estimates were available. The power TSV matched the earlier public-pair run byte for byte (SHA-256 `e140f49224fd0d59a40f3a0257ffa5bad62e385c8cf858787c576ca2c0f3357b`). |
| Report and figure | All report links resolved. The PDF contained the 79 target labels, the PNG decoded successfully, and the 79 rows in the Excel sheet matched the Dedud plotting table. |


Container release `container-v0.1.0` ran the complete public ENA pair and the complete laboratory pair in separate Slurm jobs with `--fresh`. The four ENA FASTQs passed their published MD5 checks. Each run produced 79 fresh alignments per library and matched the independently verified 316-row table byte for byte. Counts in the table below apply to both sources. [Container acceptance](container.md#acceptance-on-luminary)

The full laboratory FASTQ run used analysis release `v0.1.1` (`526d368`), CRISPResso2 2.3.4, six allocated CPUs and 18 GiB of requested memory. The paired analysis took 9 min 3 s. The run directory occupied 705 MiB. All 316 rows matched the archived paired table byte for byte. An independent calculation from the exported counts and denominators also reproduced both rates and every paired difference to within 1e-12.

| Library | Original read pairs | Pairs after preparation | Fresh target alignments |
| --- | ---: | ---: | ---: |
| NoEP | 10,084,551 | 5,971,199 | 79 |
| HiFi | 9,868,448 | 6,979,407 | 79 |

## Node-local installation

On Luminary, linking the environment on the shared filesystem was slow. The acceptance run installed the Linux packages on a compute node's local `/tmp` and ran the analysis in that allocation. The shared working directory held the FASTQs, reference index and final output.

After cloning the agent in step 1, run the Git commands in step 2 on the login node to fetch the upstream checkout. If using the public runs, also download the four files from step 3 on the login node. Defer reference indexing and manifest creation until the environment is active. Stage the locked archives in a shared cache and request an allocation from the agent repository root.

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
