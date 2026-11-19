---
name: absolve-seq
description: Analyze original ABSOLVE-seq pooled reporter FASTQs on Slurm. Select a matched donor and replicate, compare Maxi pBC criteria, export four indel-rate strategies, or reproduce manuscript treatment margins across the 79-target panel.
---

# ABSOLVE-seq analysis

Work from the repository root in the environment described in the README. Use the user's configured manifest, reference inputs and Slurm account.

1. For a donor/replicate request, run `python paired.py --list` and select the matching complete NoEP/HiFi pair. Submit `sbatch paired.sbatch --donor DONOR --replicate REPLICATE`. Use the requested Maxi method, with `manuscript` as the default.
2. For a manuscript-wide comparison, confirm that the four workbook extracts and 26 libraries are available. Use `sbatch --cpus-per-task=24 --mem=48G paired.sbatch --manuscript-panel`.
3. Inspect the scheduler's terminal state and the run output. For paired analyses verify 316 unique target/strategy rows in `results/paired_target_rates.tsv`. For the full panel verify 632 unique target/strategy/editor comparisons in `panel_comparison.tsv`.
4. Report rates in percent and differences in percentage points. For a pair, subtract NoEP from HiFi. For the manuscript comparison, use the fitted treatment margins and report the BH correction separately per strategy/editor. Distinguish threshold agreement from numerical agreement.
5. A request for filtering sensitivity may select `merged`, `any-pool`, `two-pools`, or `all-pools`. Those methods require the corresponding Maxi QC references. Use a distinct output directory for each run.
