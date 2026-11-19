# Original paired-end read processing

Choose manuscript mode to replay the original preprocessing with Trimmomatic 0.36 and paired-end CRISPResso2 2.2.13. The default for ordinary NoEP/HiFi pairs is the fast workflow. A full manuscript panel defaults to manuscript read processing. The Maxi annotation choice and the read processing mode are separate options.

## Configure the submission host

Use Python 3.10 or later, Java, Apptainer, Slurm and the upstream checkout and references from the main README. Supply the original Trimmomatic 0.36 JAR and a Linux CRISPResso2 2.2.13 image. The final panel fit uses the agent analysis image when `ABSOLVE_IMAGE` is set, or an available `Rscript` with brglm2. All paths must be visible on the compute nodes. [Trimmomatic 0.36 distribution](https://anaconda.org/bioconda/trimmomatic/files) · [CRISPResso2 2.2.13 image](https://depot.galaxyproject.org/singularity/crispresso2:2.2.13--py310hd6be1da_0)

```sh
export ABSOLVE_TRIMMOMATIC_JAR=/shared/software/trimmomatic-0.36.jar
export ABSOLVE_CRISPRESSO2213_IMAGE=/shared/images/crispresso-2.2.13.sif
export ABSOLVE_IMAGE=/shared/images/absolve-linux-amd64.sif
python panel_report/references.py /shared/submitted-workbooks data/submitted_model_reference --full-model
python paired.py --manuscript-panel --read-mode manuscript \
  --baseline-panel /shared/earlier-fast-run/panel_comparison.tsv \
  --output runs/manuscript-replay
```

Omit `--baseline-panel` to compare directly with the submitted workbooks. Choose a new run directory. Set site-specific scheduler defaults with variables such as `SBATCH_PARTITION` and `SBATCH_ACCOUNT`. The launcher submits a preparation array, a corresponding alignment array and a final aggregation job, and records their IDs in `jobs.json`. The two arrays use six CPUs, 18 GB memory and an eight-hour limit per library by default. Allow shared scratch for trimmed pairs and uncompressed per-target FASTQs.

For a single complete donor/replicate pair, run `python paired.py --donor D3734 --replicate T1 --read-mode manuscript --output runs/pair-replay`. Pair output contains 316 target/strategy rows. For a fast panel, submit `sbatch paired.sbatch --manuscript-panel --read-mode fast`.

## How reads are counted

Both input FASTQs are checked against their manifest SHA-256 hashes. Trimmomatic uses the upstream adapter file and paired clipping and quality settings. Reads are assigned using the published 23-base constant with up to five no-calls and a recognized target barcode. The 16-base pBC is recorded before removing it from R1 and from the end of R2.

For each sample and target, the worker passes both mates to CRISPResso2 2.2.13, which uses FLASH to merge them. Target batches share an amplicon, guide and fixed merge parameters across their families. Every read retains its original pBC. The worker reconstructs family allele abundances using the aligned sequences and cut-site calls. Representative-allele ties follow the v2.2.13 allele table order. Dedud applies the pinned Maxi annotation and cross-target identity check to edited alleles. Per-target counts are checked against CRISPResso’s aligned and modified totals. [Original batch command](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/absolveseq/2_demultiplex_by_plasmidBarcode.py#L89-L100)

## Inspect the results

Panel output includes 8,216 target/strategy/library count rows and 632 modeled contrasts in `results/`. The model uses the same donor-adjusted treatment margins, one-sided test and within-strategy/editor BH correction as the fast panel. With a baseline, `results/call_changes.tsv` labels each contrast as concordant in both runs, resolved, persistent or new. Agreement uses strict thresholds of Δ greater than 0.1 percentage points and q below 0.05.

Each completed library stage and target writes its output atomically. A failed alignment task can be retried for its original array index with the same run configuration. Consult the stage logs and submit the final aggregation after all affected indices complete.

## Luminary validation

The replay used 52 checksum-verified FASTQs, 26 original libraries and 79 targets. An independent check read Supplementary Tables 2 through 5, recalculated both sets of BH q values from their one-sided p values, and compared each of the 632 calls. The replay resolved 21 of the 24 fast-mode discrepancies. Three persisted and three appeared in previously concordant contrasts. Four D3734/T1 families at two targets also matched a separate per-family CRISPResso2 replay in both retained indel reads and retained total reads.

The six remaining differences comprised two Unfiltered WT calls, one Unfiltered HiFi call, two Dedup WT calls and one Dedup HiFi call. In three contrasts only the effect-size threshold changed, and in three both thresholds changed. The local comparison table records the target, editor, rates and q values for each contrast.

All 316 calls across Dedud and Dedud+Dedup agreed at the joint thresholds. Individual fitted rates and q values still differ. The replay changes trimming, read admission, paired merging and the CRISPResso2 version together, so this comparison measures their combined effect.
