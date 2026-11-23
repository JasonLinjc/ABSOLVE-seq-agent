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

## D3734/T1 Docker acceptance

The image built from `b7964ed` processed all four public D3734/T1 FASTQs from NoEP (`SRR39887646`) and HiFi (`SRR39887695`) through to `results/report.md`. The command below ran inside the new image after the input checksums were verified. No mode selection was needed. The run recorded the default paired-end CRISPResso2 2.2.13 workflow and the published plasmid barcode annotation. [Docker acceptance run](https://github.com/JasonLinjc/ABSOLVE-seq-agent/actions/runs/45116053493)

```bash
python paired.py --executor local --threads 4 --donor D3734 --replicate T1 --fresh --output /output/pair
```

We compared the outputs directly with the previous Docker run from `36b4874`, which used the same raw FASTQs. [Previous Docker acceptance run](https://github.com/JasonLinjc/ABSOLVE-seq-agent/actions/runs/45110435029)

| Output | Agreement with the previous run |
| --- | --- |
| Indel rates and paired differences | All 316 rows (79 targets × four strategies) matched byte for byte, including both arm rates and HiFi minus NoEP. All 632 underlying count rows also matched byte for byte. |
| Detection power | All 632 estimates (two assumed effects per target and strategy) matched byte for byte. |
| Figure | All 79 plotted data rows and the PNG file matched byte for byte. The PDF rendered identically and contained the same target labels. |
| Excel source data | Every cell in the `Paired Dedud` worksheet matched across all seven columns and 79 target rows. File-level differences in the PDF and Excel workbook were confined to timestamps. |
| Report | `report.md` matched byte for byte. All 11 links resolved in the downloaded report bundle. |

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
