# Run D3734/T1 from paired FASTQs

This walkthrough runs the original D3734/T1 NoEP (LV_S3) and HiFi (LV_S1) libraries through all 79 targets. Each output row includes both indel rates and their paired difference (HiFi minus NoEP). The agent trims both mates with Trimmomatic 0.36, aligns them with CRISPResso2 2.2.13 and uses the published pBC annotation and cross-target filter.

## 1. Get the agent and create an environment

Use a shared checkout and environment on a Slurm cluster. For a workstation, use the local executor in step 4. These commands install the locked Linux x86-64 environment, including Java, Trimmomatic, CRISPResso2, FLASH and the reporting tools.

```bash
git clone https://github.com/JasonLinjc/ABSOLVE-seq-agent.git
cd ABSOLVE-seq-agent
unset ABSOLVE_DATA ABSOLVE_MANIFEST ABSOLVE_CACHE ABSOLVE_RUNS ABSOLVE_IMAGE ABSOLVE_TRIMMOMATIC_JAR ABSOLVE_CRISPRESSO2213_IMAGE ABSOLVE_CRISPRESSO2213_EXECUTABLE PYTHONPATH
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source "$(conda info --base)/etc/profile.d/conda.sh"
env_prefix="$PWD/../absolve-paired-env"
mamba create -y -p "$env_prefix" --file docs/paired-linux-64.lock
conda activate "$env_prefix"
CRISPResso --version
```

Check for version 2.2.13. The launcher discovers the Trimmomatic 0.36 JAR in this environment. For platforms with compatible native packages, the dependency specification is `environment.yml`. For Linux containers or an existing cluster installation, see [container setup](container.md) and [external runtime paths](manuscript-reads.md#configure-the-runtime).

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

From the active paired environment, use a new run directory. Both mates of every selected library are checked against the manifest before trimming. The ordinary command uses the published pBC annotation.

```bash
run_dir="$PWD/runs/D3734-T1"
python paired.py --executor local --donor D3734 --replicate T1 --threads 4 --output "$run_dir"
```

### Slurm submission

With the environment and inputs visible to the compute nodes, run the following on the **submission host**. Set your site’s partition or account as needed.

```bash
export SBATCH_PARTITION=cpu1,cpu2  # Luminary example
run_dir="$PWD/runs/D3734-T1"
python paired.py --donor D3734 --replicate T1 --output "$run_dir"
cat "$run_dir/jobs.json"
```

The launcher submits preparation, alignment and reporting jobs with dependencies. Wait for the **finish job** to reach `COMPLETED` with exit code `0:0`; use its ID from `jobs.json` with `sacct -j JOB_ID --format=JobID,State,ExitCode`. Then continue with step 5.

## 5. Read the report and verify the result

Open `$run_dir/results/report.md` for target-level indel rates, HiFi minus NoEP differences and estimated detection power. The PDF, PNG and Excel source data are under `results/figures/`.

```bash
python verify_results.py pair-manuscript "$run_dir/results"
python verify_results.py power "$run_dir/results"
```

The first check compares all 316 paired results and 632 arm-level count rows against the accepted public pair. The second checks the 632 power estimates (two assumed effects for each target and strategy), the plotted values and the Excel source data. [Power method](power.md)

### Verified on Luminary

The public NoEP and HiFi FASTQs passed the ENA MD5 checks. Starting with a new reference index, the default CLI processed all 79 targets in both libraries using paired-end CRISPResso2 2.2.13 and the published annotation (finish job `21343015`). All 632 arm/target/strategy count rows matched the independent published-annotation panel recount, including numerators, denominators and retained family counts. Recomputing both arm rates and HiFi minus NoEP from those counts reproduced all 316 comparisons within 1e-12 percentage points. The power and Excel checks passed.

A second run used the local executor inside a Luminary allocation (job `21343242`). It installed the 157-package Linux lock, discovered Trimmomatic and CRISPResso2 in that environment and processed the four raw FASTQs again. The counts and indel-rate tables matched the Slurm exports byte for byte. Both runs produced the same plotted data and the same Excel values for all 79 targets.

For identical detection-power estimates, use the same locked reporting environment. With the older reporting image (NumPy 2.5.5), 11 of 632 simulated estimates differed from the lock (NumPy 1.26.4), by at most 0.9 percentage points. Recomputing power from the Slurm rate table with the lock (job `21343314`) reproduced the local executor's entire power table byte for byte.

| Accepted export | SHA-256 |
| --- | --- |
| `paired_target_rates.tsv` | `003cf3fa670477f811f5c697d0abb896525825a7991e263c1e985a78d9629a3b` |
| `library_counts.tsv` | `008a9aad4f0749d0b38c27205ca32ed2e5d58ce3aa67e73e9e56cc29d327ae9a` |
| `power_estimates.tsv` (locked environment) | `3cb3b127f0f1be573635beff772ec00e37a19361c3707fb8a341c218e9cb44e3` |

At `1617_OT_0040_ALT`, Dedud+Dedup retained 11 indel-positive HiFi families out of 691 and zero out of 1,694 for NoEP, a difference of +1.59 percentage points.

## Node-local installation

On Luminary, linking the environment on the shared filesystem was slow. The acceptance run installed the Linux packages in a compute allocation's scratch directory and ran the analysis in that allocation. The shared working directory held the FASTQs, reference index and final output.

After cloning the agent in step 1, run the Git commands in step 2 on the login node to fetch the upstream checkout. If using the public runs, also download the four files from step 3 on the login node. Defer reference indexing and manifest creation until the environment is active. Stage the locked archives in a shared cache and request an allocation from the agent repository root.

```bash
export CONDA_PKGS_DIRS="$PWD/data/conda-pkgs"
mkdir -p "$CONDA_PKGS_DIRS"
mamba create --no-rc --download-only -y -p "$PWD/data/download-env" \
    --file docs/paired-linux-64.lock
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
for entry in Path("docs/paired-linux-64.lock").read_text().splitlines():
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
python paired.py --executor local --threads 6 --donor D3734 --replicate T1 --output "$run_dir"
```
