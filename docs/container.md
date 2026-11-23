# Run the analysis container

The Linux x86-64 Dockerfile builds the paired-end analysis environment and includes the agent source, CRISPResso2 2.2.13, Python and R with brglm2. Keep FASTQs, the sample manifest, reference tables and results on mounted storage.

## Run the default paired workflow in Docker

Build from this checkout for the unified 2.2.13 runtime. The image uses the same package lock as the [paired quick start](d3734-t1-quick-start.md). With the example references and manifest prepared, mount the checkout at the same absolute path used by the manifest. Mount any external input directories at their original paths too.

```bash
docker build --platform linux/amd64 -t absolve-seq-agent:local .
run_agent() {
    docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
        --volume "$PWD:$PWD" --workdir "$PWD" absolve-seq-agent:local "$@"
}
run_agent python paired.py --executor local --donor D3734 --replicate T1 --threads 4 --output runs/D3734-T1
run_agent python verify_results.py pair-manuscript runs/D3734-T1/results
run_agent python verify_results.py power runs/D3734-T1/results
```

For an existing Linux cluster installation, the default workflow also supports Trimmomatic 0.36 on the host and a separate CRISPResso2 2.2.13 Apptainer image. [Cluster runtime configuration](manuscript-reads.md#configure-the-runtime)

## Legacy container acceptance on Luminary

Release `container-v0.1.0` used Apptainer 1.4.2 and CRISPResso2 2.3.4. The reference index built inside the SIF matched the Conda-built index file by file. The two jobs below each used six allocated CPUs and `--fresh`. Their elapsed times cover the paired analysis and exclude data transfer and indexing.

| D3734/T1 inputs | Slurm job | Elapsed | Target/strategy rows | Rows differing from the Conda baseline |
| --- | --- | --- | ---: | ---: |
| Original laboratory FASTQs | `21277384` | 9 min 14 s | 316 | 0 |
| Full ENA runs `SRR39887646` and `SRR39887695` | `21277504` | 9 min 6 s | 316 | 0 |

Both jobs processed 10,084,551 NoEP read pairs and 9,868,448 HiFi read pairs. All four strategy counts, denominators, indel rates and paired differences matched the baseline. Both exported tables had SHA-256 `0b0465ecc5aa244cc2a10654e153dbe8467322cbedd4f5bfc72806058018980b`. The verified SIF had SHA-256 `bdf50ad7f94c562b99c1bdcd869bf5b7f4c362d8f4d9590088439df52e7fb972`.

## Stage public inputs from a networked machine

When the cluster has slow external access, use the same agent revision on a networked workstation and on the cluster. Follow the [paired quick start](d3734-t1-quick-start.md) to fetch the pinned upstream checkout and four public FASTQs. Transfer the inputs to a shared cluster checkout.

```bash
destination=CLUSTER:/shared/path/to/ABSOLVE-seq-agent/
rsync -aP --relative data/upstream/test/data/ data/sra/SRP722226/ "$destination"
```

On the cluster, build the reference index and create the manifest with the `--ena` option to verify the transferred files and record their cluster paths. Then follow the [Slurm submission and verification steps](d3734-t1-quick-start.md#slurm-submission).
