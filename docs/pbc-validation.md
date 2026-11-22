# Published pBC annotation validation across 26 libraries

On Luminary, the default published pBC annotation and the annotation rebuilt from the original six Maxi plasmid pool FASTQs with `--published-maxi` produced byte-identical count and model tables **within each read mode**. The historical two-tBC exception remains a separate option (`--historical-maxi`). [Build commands and rule definitions](historical-maxi.md#build-the-reference-once)

| Read mode | Alignments shared by the two annotation runs | Library counts | Model contrasts |
| --- | --- | --- | --- |
| Fast | Archived fast-mode CRISPResso2 alignments | 8,216/8,216 identical | 632/632 identical under each of two donor-weighting schemes |
| Manuscript reproduction | One paired-end CRISPResso2 2.2.13 replay with FLASH | 8,216/8,216 identical | 632/632 identical with equal donor weights |

The 26 libraries comprised nine NoEP, nine WT and eight HiFi libraries. Each annotation run independently counted 79 targets under four strategies. The count comparison covered indel numerators, denominators and pBC family counts. Each model compared WT and HiFi with NoEP for every target and strategy. Control and edited indel rates, Δindel rates, one-sided p values, BH-adjusted q values and positive calls all matched exactly between annotations.

## Run the full acceptance check

From the repository root on Luminary, run:

```bash
sbatch --wait -p cpu1,cpu2 validate_pbc_panel.sbatch
```

The job reads `data/pbc-validation.json`, allocates 32 CPUs and writes a new `runs/pbc-validation-*/` directory each time. The first line of the Slurm log gives its absolute path. To use another configuration or a specific new output directory, append `--config /path/to/config.json --output /path/to/new-run` to the same command.

### Configure the archives once

Copy the [example configuration](../examples/pbc-validation.json) to `data/pbc-validation.json` and set the paths below. Relative paths resolve from the config file's directory. The compute node needs Python 3.9 or later and Apptainer with read access to the archives.

| Config key | Input |
| --- | --- |
| `data` | Panel data root containing `references/target_info/`, `references/target_guide_mapping.json` and the indexed default annotation in `references/published/`. |
| `fast_inventory` | Full fast-panel `sample_inventory.tsv` with the original 26 library identities, R1/R2 hashes and absolute `folder` paths. Each folder contains the prepared reads, job metadata and retained per-target CRISPResso outputs. Keep the original paths referenced by the metadata accessible. |
| `manuscript_alignments` | Shared replay directory containing `<sample_name>/<target>/alignment.json`, `CRISPResso_output.fastq.gz` and the CRISPResso quantification TSV. The receipt identifies the sample, target, input pair count and aligned-file SHA-256. |
| `manuscript_archive` | Original manuscript-mode default results directory containing `sample_inventory.tsv`, `library_counts.tsv` and `panel_comparison.tsv`. |
| `rebuilt_reference` | Saved annotation built with `setup_references.py --published-maxi` from the original six Maxi FASTQs. |
| `submitted_tables` | Submitted workbook extracts used by the panel model (`SupplementaryTable2.csv` through `SupplementaryTable5.csv`). |
| `fast_stats_image`, `manuscript_stats_image` | Apptainer images with R and brglm2 used in the respective archived fits. Matching the archived runtime matters for byte-level model checks. |

### Read the results

The command reuses the alignments, recounts both annotations independently and refits the models. It compares the two annotation exports within each read mode, checks fast exports against the accepted digests below and compares the manuscript default exports directly with the original archive.

| Output | Contents |
| --- | --- |
| `summary.tsv` | One row per comparison with the row counts, number of differences and both file digests. |
| `differences.tsv` | One row per differing field with the target, strategy, library or editor, and both values. File-level differences include row ordering or formatting. A header-only file means no differences. |
| `SHA256SUMS` | SHA-256 for the count and model exports, manuscript archive copies and both comparison tables. Run `sha256sum -c SHA256SUMS` from the results directory to verify the saved files. |
| `run.json` | Resolved input paths, worker count and run status. The per-mode directories contain the tables used in the comparisons. |

A successful check exits **0**. A completed comparison with differences exits **1**. An input or execution error exits **2** and is reported in the Slurm log. Run the command again to create a separate result directory.

## Manuscript archive check

The original manuscript-mode run retained target-split R1/R2 FASTQs and its default-annotation result tables. Its temporary per-read alignments had been cleared. The replay used the archived FASTQs, the original CRISPResso2 image and the same FLASH parameters. Both recounts read the same 2,054 library/target alignment files and verified their SHA-256 digests. Each model fit weighted the three donors equally on the logit scale.

Both replayed exports (`library_counts.tsv` and `panel_comparison.tsv`) also matched **the original manuscript-mode default archive byte for byte**, including the counts before cross-target filtering. Agreement with the submitted supplementary tables is documented separately. [Supplementary-table comparison](manuscript-reads.md#luminary-validation)

## Export checksums

Each digest below applies to the default and rebuilt annotation exports for that row. The manuscript-mode digests also apply to the archived default exports.

| Read mode | Export | SHA-256 |
| --- | --- | --- |
| Fast | Library counts | `3db1f8203892c29edd231e2d88d942223b9c6a5b44a6d6d23aabc107f90a7a84` |
| Fast | Model contrasts, equal donor weights | `627df9a661d2160184753c39a6bcf39dcef40449c42df7f867b229b1f0342755` |
| Fast | Model contrasts, pBC family weights | `ff94527fbdf0a95bcfc59da57b708ed838692368c34ca69c046c01c537bfaa99` |
| Manuscript | Library counts | `aefad456069a4eb3d8eef4e2cdef9b913c6e8d2791df2be04ea306f60b9fef40` |
| Manuscript | Model contrasts, equal donor weights | `3f52e3068390b9bf5294a816180d2dde7107a9cd29458f0f50b0079341c6290d` |
