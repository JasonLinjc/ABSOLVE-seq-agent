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

Panel output includes 8,216 target/strategy/library count rows and 632 modeled contrasts in `results/`. For each treatment, manuscript mode averages the fitted logits for D3734, D3742 and D3908 with weight 1/3 each, then applies the inverse logit. It subtracts the NoEP rate from the edited treatment rate. In `donor_weights.tsv`, `marginal_weight` is the weight used for prediction and `family_row_weight` is the observed share of pBC-family rows. The fast panel uses the observed shares. Both modes use a one-sided treatment-coefficient Wald test and BH correction within each 79-target strategy/editor group. A hit requires Δ greater than 0.1 percentage points and q below 0.05. With a fast baseline, `results/call_changes.tsv` labels contrasts as concordant, resolved, persistent or new.

Each completed library stage and target writes its output atomically. A failed alignment task can be retried for its original array index with the same run configuration. Consult the stage logs and submit the final aggregation after all affected indices complete.

### Refit an existing manuscript panel from counts

Run the following in the repository root with R and brglm2 available. Use a new destination directory and the workbook extracts from the original run. This recomputes all four strategies, both editors, treatment margins, and p and q values from the saved counts.

```sh
mkdir runs/equal-donor-refit
cp runs/manuscript-replay/results/library_counts.tsv runs/equal-donor-refit/
Rscript panel_report/manuscript_model.R runs/equal-donor-refit data/submitted_model_reference equal-donor
```

## Luminary validation

The paired-end replay used 52 checksum-verified FASTQs, 26 original libraries and 79 targets. Four D3734/T1 families at two targets matched a separate per-family CRISPResso2 replay in both retained indel reads and retained total reads. The initial family-weighted panel agreed with Supplementary Tables 2 through 5 on 626 of 632 calls.

The equal-donor refit reused the same 8,216 count rows. An independent check reconstructed all 632 treatment contrasts from the previous donor-specific logits and checked the submitted rates and BH q values against the four workbooks. All one-sided p and BH q values matched the previous paired-end panel exactly. The fast-mode regression check reproduced its archived panel and model-field tables byte for byte.

One contrast, Unfiltered HiFi at `1617_OT_0356_ALT`, moved from Δ 0.099694 to 0.103867 percentage points (submitted 0.106626) with q 5.10 × 10⁻¹⁰. This restored the submitted hit. No previously concordant call became discordant. The five remaining differences are listed below. Rates are rounded for display, and calls use full precision.

| Strategy | Editor | Target | Equal-donor Δ (pp) | Submitted Δ (pp) | Hit in |
| --- | --- | --- | ---: | ---: | --- |
| Unfiltered | WT | `1450_OT_0163_ALT` | 0.088032 | 0.110287 | Submitted |
| Unfiltered | WT | `1450_OT_0176_ALT` | 0.088096 | 0.130761 | Submitted |
| Dedup | WT | `1450_OT_0142_ALT` | 0.060436 | 0.110040 | Submitted |
| Dedup | WT | `1450_OT_0210_ALT` | 0.061957 | 0.118504 | Submitted |
| Dedup | HiFi | `1617_OT_0041_REF` | 0.119431 | 0.085360 | Equal-donor replay |

The final panel agreed on 627 of 632 calls, including all 316 Dedud and Dedud+Dedup contrasts. Agreement here refers to the joint Δ and q thresholds.
