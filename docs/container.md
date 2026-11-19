# Run the analysis container

The Linux x86-64 image includes the agent source, CRISPResso2, Python and R with brglm2. The Docker archive and Apptainer SIF contain the same build. FASTQs, the sample manifest, reference tables and results stay on mounted storage.

## Apptainer on Slurm

Use a shared checkout visible to the login and compute nodes. From its root, download the SIF once and verify its checksum.

```bash
release=container-v0.1.0
base="https://github.com/JasonLinjc/ABSOLVE-seq-agent/releases/download/$release"
mkdir -p data/containers
(cd data/containers && curl -fLO "$base/absolve-linux-amd64.sif" && curl -fLO "$base/SHA256SUMS" && \
    awk '$2 == "absolve-linux-amd64.sif"' SHA256SUMS | sha256sum -c -)
export ABSOLVE_IMAGE="$PWD/data/containers/absolve-linux-amd64.sif"
export ABSOLVE_DATA="$PWD/data"
export APPTAINERENV_ABSOLVE_DATA="$ABSOLVE_DATA"
apptainer exec --cleanenv "$ABSOLVE_IMAGE" CRISPResso --version
```

Use the Git and download commands in the [D3734/T1 quick start](d3734-t1-quick-start.md) for the upstream checkout and the four ENA FASTQs. If starting with a fresh reference directory and manifest, run the following commands from the shared checkout. Bind other input or output directories at the same absolute paths when they are outside the checkout. For example, set `APPTAINER_BINDPATH=/shared/fastqs:/shared/fastqs:ro` for external reads.

```bash
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python /opt/absolve/setup_references.py
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python /opt/absolve/examples/d3734_t1_manifest.py \
    --ena "$PWD/data/sra/SRP722226" "$PWD/data/sample_manifest.tsv"
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python /opt/absolve/paired.py --list
sbatch --partition=cpu1,cpu2 paired.sbatch --donor D3734 --replicate T1 --fresh \
    --output "$PWD/runs/container-D3734-T1"
```

Set the partition for your cluster. The Slurm wrapper passes the allocation and analysis settings to Apptainer and executes the code shipped inside the SIF. Use the result checks in [step 5](d3734-t1-quick-start.md#5-verify-the-result). The same wrapper accepts `--manuscript-panel` with the [full panel inputs](../README.md#manuscript-comparison).

## Docker

Download `absolve-linux-amd64.tar.gz` and `SHA256SUMS` from the [container release](https://github.com/JasonLinjc/ABSOLVE-seq-agent/releases/tag/container-v0.1.0). Verify the archive and load it. Use the exact image tag printed by `docker load` in subsequent commands.

```bash
awk '$2 == "absolve-linux-amd64.tar.gz"' SHA256SUMS | sha256sum -c -
docker load -i absolve-linux-amd64.tar.gz
```

To build locally from the repository root and check the available CLI options:

```bash
docker build --platform linux/amd64 -t absolve-seq-agent:local .
docker run --rm --platform linux/amd64 absolve-seq-agent:local
```
