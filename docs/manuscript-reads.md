# Original paired-end read processing

Both matched NoEP/HiFi pairs and the full study panel use Trimmomatic 0.36, paired-end CRISPResso2 2.2.13 and the published pBC annotation by default. To reproduce the submitted workbooks, explicitly select `--manuscript-panel --read-mode historical`. For routine analysis, you can also select a saved Maxi annotation independently of the read mode.

## Analyze all study samples

The full panel covers 26 libraries from D3734, D3742 and D3908. The design has nine NoEP, nine WT and eight HiFi libraries. D3734/T3 includes WT and NoEP only. All 26 contribute to the donor-adjusted model.

Run these Bash commands from the current agent checkout on a machine with internet access. Use the [locked Linux Conda environment](d3734-t1-quick-start.md#1-get-the-agent-and-create-an-environment) or the [Docker or Apptainer setup](container.md#download-a-prebuilt-image). These environments include R and brglm2 for the model fit. With a container, prefix each Python command below with `run_agent`.

**1. Build the published reference index.** From the checkout, fetch the five files the analysis uses from the pinned upstream revision. For a new run the commands below use the default `data/` directory.

```bash
set -e
upstream=https://raw.githubusercontent.com/pellinlab/ABSOLVE-seq/271355dfe40cf46b952e82b88dd4521f23f9feb4/test/data
ref=data/upstream/test/data
mkdir -p "$ref/target_info"
for file in adapters.fa plasmid_barcode_category.tsv.gz target_info/NovaSeq3_sample_info.csv target_info/OT_guide_amplicon_seq.csv target_info/target_oligos_sequenes.csv; do
  curl -fLsS --retry 4 --retry-all-errors "$upstream/$file" -o "$ref/$file"
done
python setup_references.py
```

**2. Prepare the 52 FASTQs.** The public FASTQs total 51.7 GB. A full run retained about 111 GiB of intermediate files in addition to the inputs. Allow extra free space for the environment and temporary alignment files. The download table provides both mate URLs and ENA MD5 checksums for each library. Use a networked workstation to stage the files if your cluster restricts internet access. [Study download table](../examples/study_panel_ena.tsv) · [Public study](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP722226)

```bash
fastqs="$PWD/data/sra/study-panel"
mkdir -p "$fastqs"
awk -F '\t' 'NR > 1 {print $4; print $6}' examples/study_panel_ena.tsv |
  (cd "$fastqs" && xargs -P 4 -n 1 curl -fLsS --retry 4 --retry-all-errors -C - -O)
python examples/d3734_t1_manifest.py --panel --ena "$fastqs" data/sample_manifest.tsv
```

The manifest builder checks each compressed FASTQ against its ENA MD5 and records a local SHA-256. Run it where the analysis will execute so the absolute paths resolve. If you already have the original laboratory FASTQs under `WT/`, `NoEP/` and `HiFi/`, use `python examples/d3734_t1_manifest.py --panel /path/to/libraries data/sample_manifest.tsv`. An existing full manifest can also be selected with `ABSOLVE_MANIFEST` in a native environment.

**3. Run the panel and verify the report.** On a Slurm submission host, load the scheduler commands and confirm `command -v sbatch` succeeds. The standard workflow fits all editing effects and estimates power directly from the sequencing counts. For an optional comparison with the paper, place the four full-model reference extracts and their `source_sha256.json` in `data/submitted_model_reference/` before launch.

```bash
python paired.py --manuscript-panel --output runs/study-panel
```

To use a workstation or one existing allocation, add `--executor local --threads N` with the allocated CPU count. Local execution returns after the report is written. On Slurm, wait for the finish job recorded in `runs/study-panel/jobs.json`. Then run `python verify_results.py power runs/study-panel/results` and open `runs/study-panel/results/report.md`. [Result files](../README.md#study-panel-results)

### New-user acceptance

We checked the public run IDs, mate URLs and archive MD5 values for all 26 libraries against ENA, tested all 52 download links and verified the four D3734/T1 files in full. The full-panel analysis used staged original laboratory FASTQs. Each library's input read-pair count matched its ENA record, totaling 334,970,765 pairs. We also downloaded all five pinned public reference files and rebuilt the 79-target index.

On Luminary we built a fresh environment, index and sample manifest, then ran the routine workflow from FASTQs with no supplementary tables or cached alignments. All 2,054 library/target alignments completed. The report contained 8,216 count rows, 632 WT or HiFi contrasts, 1,264 power estimates and 158 Figure 3b source rows. The output verifier checked the plotted values and Excel cells. All 14 report links resolved.

Against the earlier full-panel run with the supplementary tables, the library counts were byte-identical (SHA-256 `aefad456069a4eb3d8eef4e2cdef9b913c6e8d2791df2be04ea306f60b9fef40`). All 632 contrasts matched for indel rates, editing differences, p-values, q-values and calls. All 18,012 fitted model fields matched as well.

For exact simulation results across runs, use the same runtime. Within each tested runtime, adding the supplementary tables produced identical power tables. Between the locked Conda environment (NumPy 1.26.4) and the v0.2.0 release container (NumPy 2.5.5), 79 of 1,264 simulated power estimates differed, by at most 1.2 percentage points. One power cell in the Figure 3b source data changed from 32.6% in the container to 32.5% in Conda. The container reproduced the previous run's power table and plotted data byte for byte, and every Figure 3b data cell matched.

## Reproduce the submitted calls from the original FASTQs

> I have all 26 libraries, their manifest and submitted Supplementary Tables 2–5. Check the historical software versions, reproduce the submitted tables and summarize agreement for WT and HiFi by strategy. Estimate power from my counts, draw Figure 3b from the fitted Dedud estimates and simulated power, and return the report, PDF and Excel source data.

Select `--manuscript-panel --read-mode historical --fresh` to process all 52 raw FASTQs. Each library goes through input SHA-256 verification, trimming, target assignment, paired-end alignment and counting in a new run directory. The panel uses the historical six-Maxi pBC membership and the existing cross-target filter.

| Libraries | Target assignment after historical trimming | Paired reads and tied family representatives |
| --- | --- | --- |
| 8 HiFi and 9 NoEP | Exact tBC at offsets 39–47 and first occurrence of the 16-base constant motif at offset 16 in R1 (zero-based) | CRISPResso2 2.2.13 with single-threaded FLASH. Among equally abundant alleles, use the first encountered in the aligned FASTQ, including after cross-target filtering. |
| 9 WT | Same tBC position and the 23-base constant rule (up to five no-calls) | CRISPResso2 2.3.2 with its default fastp 0.24.1 merger. Break abundance ties by aligned read and then aligned reference. |

Configure the [runtime](#configure-the-runtime), the 26-library manifest and the four full-model workbook extracts. Set two additional paths. The WT environment must contain its working `bin/python3.12`, `bin/CRISPResso` and `bin/fastp`. The pBC directory comes from the [historical six-Maxi build](historical-maxi.md#build-the-reference-once).

```bash
export ABSOLVE_WT232_ENV=/shared/envs/crispresso2-2.3.2
export ABSOLVE_HISTORICAL_MAXI_REFERENCE=/shared/references/historical-maxi
python paired.py --manuscript-panel --read-mode historical --fresh --output runs/historical-panel
```

On Slurm, follow the three job IDs in `runs/historical-panel/jobs.json` through the finish job. For a local run add `--executor local --threads N`. When the report is ready, compare the calls directly against the four original workbooks and verify the power table, Figure 3b and its Excel source.

```bash
python verify_results.py panel-historical runs/historical-panel/results --workbooks /path/to/submitted/workbooks
python verify_results.py power runs/historical-panel/results
```

The workbook check independently computes BH-adjusted q-values within each editor and strategy. It checks all 632 calls and reports the largest absolute differences in the control rate, edited rate, Δ and q. The acceptance tolerances are 0.01 percentage points for rates and Δ, and 0.01 for q. The per-comparison values remain in `panel_comparison.tsv`.

| Starting point | Effect estimate | Report path relative to `--output` |
| --- | --- | --- |
| Matched pair, uniform paired-end reads | Observed HiFi and NoEP rates for one donor and replicate | `results/report.md` |
| Full panel, `--read-mode historical --fresh` | Arm-specific historical processing from the 52 FASTQs and equal donor weights | `results/report.md` |
| Full panel, explicit `--read-mode manuscript` | Three donors weighted equally on the logit scale | `results/report.md` |
| Full panel, `--refit-counts FILE --output DIR` | Three donors weighted equally on the logit scale, using the filtering already present in the counts | `report.md` |

### Full-panel acceptance on Luminary

[Results against Supplementary Tables 2–5, the residual at `1450_OT_0235_ALT`, and export checksums](historical-panel-acceptance.md)

## Collect a historical per-pBC batch

Use this entry point when the starting material is one archived CRISPRessoBatch directory with a separate CRISPResso run per plasmid barcode. Supply the original library ID and group label from the sample metadata, and choose a new output directory. This example uses the WT D3908/T1 batch. The collector requires pandas (tested with 2.2.3).

```bash
python -m panel_report.collect_historical_families \
  --batch /path/to/CRISPRessoBatch_on_WT_D3908_T1_1450-OT-0235-ALT_input \
  --library-id LV-S24 --group WTCas9 --output runs/D3908-T1-WT-0235/alleles
```

The collector compares the unique pBCs in the `UMI` column of each combined table with `results.completed_batch_arr` in `CRISPResso2Batch_info.json`. It checks the guide-window and full-length tables independently, before computing read percentages. Exit 0 publishes both historical `*withUMI.txt` tables together for downstream counting. A membership mismatch exits 1, lists each affected filename and its missing or unexpected pBCs, and leaves the output directory uncreated. After restoring the inputs, rerun the same command. [Three-family restart control](historical-panel-acceptance.md#three-family-restart-control-on-luminary)

## Configure the runtime

Use Python 3.10 or later, Java and the upstream checkout and references from the main README. Use the native paired environment either locally or on a Slurm cluster. Linux hosts with Apptainer can also run the alignment in a separate image. The paired Conda environment includes Trimmomatic 0.36, FLASH and native CRISPResso2 2.2.13. The launcher finds its JAR in the environment. To use an external JAR and Linux alignment image, set the paths below. The final panel fit uses the agent analysis image when `ABSOLVE_IMAGE` is set, or an available `Rscript` with brglm2. All paths must be visible on the compute nodes. [Trimmomatic 0.36 distribution](https://anaconda.org/bioconda/trimmomatic/files) · [CRISPResso2 2.2.13 image](https://depot.galaxyproject.org/singularity/crispresso2:2.2.13--py310hd6be1da_0)

```sh
export ABSOLVE_TRIMMOMATIC_JAR=/shared/software/trimmomatic-0.36.jar
export ABSOLVE_CRISPRESSO2213_IMAGE=/shared/images/crispresso-2.2.13.sif
export ABSOLVE_IMAGE=/shared/images/absolve-linux-amd64.sif
python panel_report/references.py /shared/submitted-workbooks data/submitted_model_reference --full-model
python paired.py --manuscript-panel --read-mode manuscript \
  --baseline-panel /shared/previous-panel/panel_comparison.tsv \
  --output runs/manuscript-replay
```

Omit `--baseline-panel` to compare directly with the submitted workbooks. Choose a new run directory. On Slurm, set site-specific scheduler defaults with variables such as `SBATCH_PARTITION` and `SBATCH_ACCOUNT`. The launcher submits a preparation array, a corresponding alignment array and a final aggregation job, and records their IDs in `jobs.json`. The two arrays use six CPUs, 18 GB memory and an eight-hour limit per library by default. Allow shared scratch for trimmed pairs and uncompressed per-target FASTQs.

For a single complete donor/replicate pair, run `python paired.py --donor D3734 --replicate T1 --read-mode manuscript --output runs/pair-replay`. Pair output contains 316 target/strategy rows.

## Run locally

When neither a Slurm allocation nor `sbatch` is available, the same Python commands process each library on this machine and return after aggregation. Use `--threads N` to set the trimming and target alignment concurrency within one library. Libraries run sequentially. Allocate working space for trimmed FASTQs and per-target alignments across all requested libraries.

For the image route, keep the JAR and image settings above. For a native run, use the following environment with the original JAR, a working CRISPResso2 2.2.13 executable, and FLASH on PATH. For a full panel, also configure R with brglm2. For a submitted-workbook comparison, prepare the full-model reference extracts. The workflow checks the historical tool versions.

```sh
export ABSOLVE_TRIMMOMATIC_JAR=/software/trimmomatic-0.36.jar
unset ABSOLVE_CRISPRESSO2213_IMAGE ABSOLVE_IMAGE
export PATH="/software/crispresso-2.2.13/bin:$PATH"
python paired.py --executor local --read-mode manuscript --donor D3734 --replicate T1 \
  --threads 4 --output runs/local-pair
python paired.py --executor local --manuscript-panel --threads 4 --output runs/local-panel
```

To select a specific native CRISPResso executable, set `ABSOLVE_CRISPRESSO2213_EXECUTABLE` to its path. To analyze frozen counts directly, use [the count refit](historical-maxi.md#refit-existing-counts). Reference rebuilding and Figure 3b export also run directly from Python.

## How reads are counted

Both FASTQs are checked against the manifest SHA-256 hashes. Paired-end Trimmomatic 0.36 uses Phred+33 and the upstream adapter file with `ILLUMINACLIP:adapters.fa:2:30:10`. The remaining steps follow the processing record for each original library, in the order shown.

| Libraries | Trimmomatic steps after adapter clipping |
| --- | --- |
| HiFi (8) and NoEP (9) | `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 AVGQUAL:20 MINLEN:48 CROP:105` |
| WT (9) | `LEADING:3 TRAILING:3 AVGQUAL:20 MINLEN:48` |

Each of the 17 HiFi and NoEP input filenames matched a completed command in the archived trimming log. The WT archive contains several trial commands. For all nine WT libraries, the pair count recorded at the start of target demultiplexing matched the run without a sliding window, and the number of admitted reads matched the separate WT target-count table. The selected steps are saved per library in `run.json` and `results/sample_inventory.tsv`. Configurations created before historical trimming was added retain their original steps when resumed.

Only paired survivors enter target assignment. Reads are assigned using the published 23-base constant with up to five no-calls and a recognized target barcode. The 16-base pBC is recorded before removing it from R1 and from the end of R2.

For each sample and target, the worker passes both mates to CRISPResso2 2.2.13, which uses FLASH to merge them. Target batches share an amplicon, guide and fixed merge parameters across their families. Every read retains its original pBC. The worker reconstructs family allele abundances using the aligned sequences and cut-site calls. Representative-allele ties follow the v2.2.13 allele table order. Dedud selects pBCs from the configured reference (the published Maxi annotation by default), then applies the cross-target identity check to edited alleles. Per-target counts are checked against CRISPResso’s aligned and modified totals. [Original batch command](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/absolveseq/2_demultiplex_by_plasmidBarcode.py#L89-L100)

## Inspect the results

Panel output includes 8,216 target/strategy/library count rows and 632 modeled contrasts in `results/`. For each treatment, manuscript mode averages the fitted logits for D3734, D3742 and D3908 with weight 1/3 each, then applies the inverse logit. It subtracts the NoEP rate from the edited treatment rate. In `donor_weights.tsv`, `marginal_weight` is the weight used for prediction and `family_row_weight` is the observed share of pBC-family rows. The panel uses a one-sided treatment-coefficient Wald test and BH correction within each 79-target strategy/editor group. A hit requires Δ greater than 0.1 percentage points and q below 0.05. With a previous panel as the baseline, `results/call_changes.tsv` labels contrasts as concordant, resolved, persistent or new. Its legacy `fast_*` columns contain the values from the supplied baseline.

Each completed library stage and target writes its output atomically. A failed alignment task can be retried for its original array index with the same run configuration. Consult the stage logs and submit the final aggregation after all affected indices complete.

### Refit an existing manuscript panel from counts

Run the following in the repository root with R and brglm2 available. Use a new destination directory and the workbook extracts from the original run. This recomputes all four strategies, both editors, treatment margins, and p and q values from the saved counts. To apply a change to trimming, launch a fresh manuscript panel from the FASTQs.

```sh
mkdir runs/equal-donor-refit
cp runs/manuscript-replay/results/library_counts.tsv runs/equal-donor-refit/
Rscript panel_report/manuscript_model.R runs/equal-donor-refit data/submitted_model_reference equal-donor
```

## Luminary validation

The historical-trimming replay started from 52 checksum-verified FASTQs, realigned all 26 libraries across 79 targets, and exported 8,216 count rows. For every library, all five Trimmomatic counts matched its historical log for input pairs, paired survivors, R1-only survivors, R2-only survivors and dropped pairs. The WT target-admission totals also matched the archived downstream counts for all nine libraries.

An independent check compared the new equal-donor margins and BH q values with Supplementary Tables 2 through 5. Agreement rose from 627 to **628 of 632 calls** (99.37%). One of the five previously discordant calls resolved, four persisted and none appeared among the 627 previously concordant calls. All 316 HiFi calls agreed, as did all Dedud and Dedud+Dedup calls for both editors.

At `1617_OT_0041_REF` (Dedup, HiFi), six previously admitted singleton pBC families dropped out during historical trimming. The HiFi edited representatives fell from 55/18,825 to 49/18,229 families. The public Dedup allele tables contained 48/18,290. The two reads from the separately traced tied family still passed. The modeled Δ fell from 0.119431 to 0.092321 percentage points and q rose from 0.013139 to 0.068430. Both the replay and the submitted table (Δ 0.085360, q 0.101886) classified the site below the joint hit threshold.

All four remaining discrepancies were WT hits in the submitted tables. Values below are rounded for display.

| Strategy | Target | Replay Δ (pp) | Submitted Δ (pp) | Replay q | Submitted q |
| --- | --- | ---: | ---: | ---: | ---: |
| Unfiltered | `1450_OT_0163_ALT` | 0.088129 | 0.110287 | 0 | 0 |
| Unfiltered | `1450_OT_0176_ALT` | 0.088283 | 0.130761 | 2.43 × 10⁻⁹ | 0 |
| Dedup | `1450_OT_0142_ALT` | 0.057495 | 0.110040 | 0.444715 | 0.028716 |
| Dedup | `1450_OT_0210_ALT` | 0.061067 | 0.118504 | 0.091502 | 0.000254 |

The two Unfiltered contrasts crossed only the Δ threshold. The two Dedup contrasts crossed both thresholds. WT library count rows matched the previous replay exactly. The shared model also uses the newly trimmed NoEP and HiFi counts.
