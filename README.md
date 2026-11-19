# ABSOLVE-seq agent

A Codex skill and Slurm workflow for measuring cut-site indels in the original ABSOLVE-seq pooled reporter libraries. Choose a donor, replicate and Maxi plasmid-barcode (pBC) rule to get four indel-rate estimates for every target. The manuscript mode compares WT and HiFi against NoEP across the original 26 libraries and 79 targets.

The analysis is based on the [upstream ABSOLVE-seq repository](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4). Read the [upstream license](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/LICENSE.md) for its code and reference assets.

## Quick start

Try the agent with public ABSOLVE-seq data from [SRA study SRP722226](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP722226) (BioProject PRJNA1503690). The D3734/T1 example uses just two paired-end runs, NoEP `SRR39887646` and HiFi `SRR39887695`. Follow the [D3734/T1 quick start](docs/d3734-t1-quick-start.md) to download the four FASTQs, verify their archive checksums and export the 316-row comparison.

Clone this repository and run commands from its root. The analysis runs on a Linux Slurm cluster with a shared filesystem. Activate the environment **before** submitting jobs. Supply your site's partition, account or wall time through the usual `sbatch` options.

```sh
conda env create -f environment.yml
conda activate absolve-agent
python paired.py --help
ABSOLVE_MANIFEST=examples/manifest.example.tsv python paired.py --list
```

The example manifest demonstrates pair discovery. For analysis, provide a tab-delimited manifest at `data/sample_manifest.tsv` with absolute R1/R2 paths, their SHA-256 values and the original sample names. The quick start generates this manifest from either the public D3734/T1 files in their ENA filenames or the laboratory files in their original filenames. See the [manifest template](examples/manifest.example.tsv). The workflow checks file identity before processing. If you have the laboratory inventory and SRA submission workbook, generate the manifest with `python resolve_inputs.py INVENTORY.tsv METADATA.xlsx data/sample_manifest.tsv`.

Fetch the public sample table, target definitions and published Maxi annotation from the pinned revision. Indexing writes target-specific files into a new `data/references/` directory.

```sh
mkdir -p data
git clone https://github.com/pellinlab/ABSOLVE-seq.git data/upstream
git -C data/upstream checkout --detach 271355dfe40cf46b952e82b88dd4521f23f9feb4
python setup_references.py
python paired.py --list
sbatch paired.sbatch --donor D3734 --replicate T1
```

In Codex, open the repository and use the [absolve-seq skill](.agents/skills/absolve-seq/SKILL.md). For example, ask “Analyze D3734/T1 with the default manuscript method. Export the NoEP and HiFi indel rates and their difference for each of the four counting strategies.”

## Paired results

Each completed pair produces `runs/paired_JOBID/results/paired_target_rates.tsv`, with 79 targets × 4 strategies. Each row contains counts, denominators, percent rates and HiFi minus NoEP in percentage points. A rate with a zero denominator is blank. Use `--output DIRECTORY` for a new destination.

| Strategy | Counting unit | Filtering |
| --- | --- | --- |
| Unfiltered | Aligned reads | All aligned families |
| Dedup | One representative allele per sample–target–pBC | Most frequent allele, with top ties resolved by first occurrence |
| Dedud | Aligned reads | Good pBC and cross-target identity check |
| Dedud+Dedup | One representative per retained family | Select the representative after filtering |

The default `--maxi-method manuscript` uses the published pBC annotation. Edited alleles must match the assigned target uniquely by minimum edit distance across the panel. Use `--target-filter pbc-only` to replay the earlier barcode-only counts. The four alternative Maxi methods use locally supplied QC outputs.

| Maxi method | Evidence for a Good pBC |
| --- | --- |
| `merged` | QC of pooled Maxi reads |
| `any-pool` | Pass in at least one independent pool |
| `two-pools` | Pass in at least two independent pools |
| `all-pools` | Pass in all three independent pools |

Set `ABSOLVE_MAXI_MERGED` to the merged reference directory and `ABSOLVE_MAXI_POOLS` to the parent of `LV_Maxi_1`, `LV_Maxi_2` and `LV_Maxi_3`. Each reference requires `summary.json` and its checksum-matched `plasmid_barcode_category.tsv.gz`. For example, submit `sbatch paired.sbatch --donor D3742 --replicate T2 --maxi-method two-pools --min-exact-reads 3`.

## Manuscript comparison

Provide the submitted files `SupplementaryTable2.xlsx` through `SupplementaryTable5.xlsx` and all 26 original libraries. The sample design includes the WT and NoEP libraries for D3734/T3, which has no matched HiFi library.

```sh
python panel_report/references.py /path/to/submitted/workbooks data/submitted_model_reference --full-model
sbatch --cpus-per-task=24 --mem=48G paired.sbatch --manuscript-panel
```

The model uses Firth binomial regression with treatment and donor, the observed donor frequencies among eligible pBC-family rows for treatment margins on the logit scale, and a one-sided Wald test with BH correction across 79 targets per editor and strategy. The primary comparison uses Δ greater than 0.1 percentage points and BH q below 0.05.

| Output under `runs/manuscript_JOBID/` | Contents |
| --- | --- |
| `comparison_summary.tsv` | Threshold agreement for each strategy and editor |
| `panel_comparison.tsv` | 632 comparisons with agent and submitted rates, differences, p and q |
| `model_field_comparison.tsv` | 57 submitted model fields per target and strategy |
| `library_counts.tsv` | Per-library counts before and after cross-target filtering |
| `donor_weights.tsv`, `sample_inventory.tsv`, `provenance.json` | Marginalization weights, library inventory and source hashes |

## Local inputs and caches

`ABSOLVE_DATA` overrides the default `data/` directory. `ABSOLVE_MANIFEST` selects an individual manifest. `ABSOLVE_CACHE` can point to a compatible historical deployment containing paired or WT alignment provenance. With fresh inputs the agent prepares and aligns libraries with CRISPResso2. Allow sufficient cluster scratch and wall time for your FASTQs. On a single pair, `--fresh` forces realignment. `Rscript` must find `brglm2` through the active environment or `R_LIBS_USER`.

The upstream reference, sample manifest, FASTQs, submitted workbooks and run outputs live outside version control. The scripts in `panel_report/` also support inspecting historical matched cohorts, including pBC family tracing; the archived-cohort audit scripts expect the original provenance and alignment layout.
