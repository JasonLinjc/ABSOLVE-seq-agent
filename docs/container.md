# Run the analysis container

The Linux x86-64 image includes the agent source, CRISPResso2, Python and R with brglm2. The Docker archive and Apptainer SIF contain the same build. FASTQs, the sample manifest, reference tables and results stay on mounted storage.

## Acceptance on Luminary

Release `container-v0.1.0` used Apptainer 1.4.2 and CRISPResso2 2.3.4. The reference index built inside the SIF matched the Conda-built index file by file. The two jobs below each used six allocated CPUs and `--fresh`. Their elapsed times cover the paired analysis and exclude data transfer and indexing.

| D3734/T1 inputs | Slurm job | Elapsed | Target/strategy rows | Rows differing from the Conda baseline |
| --- | --- | --- | ---: | ---: |
| Original laboratory FASTQs | `21277384` | 9 min 14 s | 316 | 0 |
| Full ENA runs `SRR39887646` and `SRR39887695` | `21277504` | 9 min 6 s | 316 | 0 |

Both jobs processed 10,084,551 NoEP read pairs and 9,868,448 HiFi read pairs. All four strategy counts, denominators, indel rates and paired differences matched the baseline. Both exported tables had SHA-256 `0b0465ecc5aa244cc2a10654e153dbe8467322cbedd4f5bfc72806058018980b`. The verified SIF had SHA-256 `bdf50ad7f94c562b99c1bdcd869bf5b7f4c362d8f4d9590088439df52e7fb972`.

## Apptainer on Slurm

Use a shared checkout visible to the login and compute nodes. For a scripted SSH session, start a login shell (for example, `ssh CLUSTER 'bash -lc "command -v sbatch && command -v apptainer"'`). On Luminary, the login profile adds Slurm to the command path. Use a current checkout for the power and figure steps. The versioned SIF supplies the scientific environment.

```bash
git clone https://github.com/JasonLinjc/ABSOLVE-seq-agent.git
cd ABSOLVE-seq-agent
```

From the checkout root, download the SIF once and verify its checksum.

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
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python "$PWD/setup_references.py"
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python "$PWD/examples/d3734_t1_manifest.py" \
    --ena "$PWD/data/sra/SRP722226" "$PWD/data/sample_manifest.tsv"
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python "$PWD/paired.py" --list
job_id=$(sbatch --parsable --partition=cpu1,cpu2 paired.sbatch --donor D3734 --replicate T1 --fresh \
    --output "$PWD/runs/container-D3734-T1")
sacct -j "$job_id" --format=JobID,State,ExitCode,Elapsed
```

Set the partition for your cluster. At the tagged release, the Slurm wrapper executes the code shipped inside the SIF. The current wrapper runs the shared checkout with the same container runtime. Use the result checks in [step 5](d3734-t1-quick-start.md#5-read-the-report-and-verify-the-result). For the original paired-end workflow, use the current checkout and the host-side [manuscript submission command](manuscript-reads.md). The locked image supplies R for the final model fit.

Wait for the job to reach Slurm state `COMPLETED` with exit code `0:0`. In the shared checkout, run both checks with the same image as the job. This also works when the login node has no Python plotting or Excel packages.

```bash
results="$PWD/runs/container-D3734-T1/results"
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python "$PWD/verify_results.py" pair-fast "$results"
apptainer exec --cleanenv --bind "$PWD" "$ABSOLVE_IMAGE" python "$PWD/verify_results.py" power "$results"
```

### Stage public inputs from a networked machine

If the cluster has slow external access, use a networked workstation for the downloads. Use the same agent revision on both machines. On the workstation, follow the download commands above for the SIF and `SHA256SUMS`, and the [paired quick start](d3734-t1-quick-start.md) for the pinned upstream checkout and the four public FASTQs. From the workstation checkout, replace the destination below with the SSH host and absolute path of the new cluster checkout.

```bash
destination=CLUSTER:/shared/path/to/ABSOLVE-seq-agent/
rsync -aP --relative data/containers/ data/upstream/test/data/ data/sra/SRP722226/ "$destination"
```

Back in the cluster checkout, verify the SIF against `SHA256SUMS`, set the image and data environment variables above, then run the reference builder and the `--ena` manifest helper there. Building the manifest on the cluster records paths visible to the job and verifies the ENA MD5 of each transferred mate. Submit the documented `--fresh` command and verify the results in the image.

## Docker

Download `absolve-linux-amd64.tar.gz` and `SHA256SUMS` from the [container release](https://github.com/JasonLinjc/ABSOLVE-seq-agent/releases/tag/container-v0.1.0). Verify the archive, load it and check the analysis runtime.

```bash
awk '$2 == "absolve-linux-amd64.tar.gz"' SHA256SUMS | sha256sum -c -
docker load -i absolve-linux-amd64.tar.gz
image=absolve-seq-agent:2fde728f8359d79cd17eceb161e53aaced9b67f8
docker run --rm --platform linux/amd64 "$image" CRISPResso --version
```

To build locally, run these commands from the repository root. The second command displays the available CLI options.

```bash
docker build --platform linux/amd64 -t absolve-seq-agent:local .
docker run --rm --platform linux/amd64 absolve-seq-agent:local
```

From a current checkout, use the loaded image as a local runtime for the fast workflow. Prepare the inputs and manifest following the [paired quick start](d3734-t1-quick-start.md). Mount the checkout at its original absolute path so the manifest and the output refer to the same files on the host and inside the container. Mount any external FASTQ directories at their original paths as well.

```bash
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
    --volume "$PWD:$PWD" --workdir "$PWD" "$image" \
    python paired.py --donor D3734 --replicate T1 --threads 4 --fresh --output runs/local-pair
python verify_results.py pair-fast runs/local-pair/results
```

For a locally built image use `image=absolve-seq-agent:local`. The [manuscript workflow](manuscript-reads.md#run-locally) uses its own historical read-processing runtime.
