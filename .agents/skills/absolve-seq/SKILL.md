---
name: absolve-seq
description: Analyze original ABSOLVE-seq pooled reporter FASTQs on Slurm. Select a matched donor and replicate, compare Maxi pBC criteria, export four indel-rate strategies, or reproduce manuscript treatment margins across the 79-target panel.
---

# ABSOLVE-seq analysis

Work from the repository root in the environment described in the README. Use the user's configured manifest, reference inputs and Slurm account.

For a public test dataset, follow the [D3734/T1 quick start](../../../docs/d3734-t1-quick-start.md#public-test-data-from-sra). Download both mates of NoEP `SRR39887646` and HiFi `SRR39887695` from study `SRP722226`. Use the example's `--ena` manifest option to verify the four files before submitting the pair with `--fresh`.

1. For a donor/replicate request, run `python paired.py --list` and select the matching complete NoEP/HiFi pair. Submit `sbatch paired.sbatch --donor DONOR --replicate REPLICATE`. Use the requested Maxi method, with `manuscript` as the default.
2. For a manuscript-wide comparison, confirm that the four workbook extracts and 26 libraries are available. Follow the [paired-end reproduction guide](../../../docs/manuscript-reads.md) and run `python paired.py --manuscript-panel --read-mode manuscript` from a Slurm submission host. The run records the historical trimming steps for each library. To track changes against an earlier panel, supply its `panel_comparison.tsv` with `--baseline-panel`.
3. Inspect the scheduler's terminal state and the run output. For paired analyses verify 316 unique target/strategy rows in `results/paired_target_rates.tsv`. For the full panel verify 632 unique target/strategy/editor comparisons in the reported results directory’s `panel_comparison.tsv`.
4. Report rates in percent and differences in percentage points. For a pair, subtract NoEP from HiFi. For the manuscript panel, average the three donor-specific logits equally for each treatment before converting to rates. Report BH correction separately per strategy/editor. Distinguish threshold agreement from numerical agreement.
5. A request for filtering sensitivity may select `merged`, `any-pool`, `two-pools`, or `all-pools`. Those methods require the corresponding Maxi QC references. Use a distinct output directory for each run.

To rebuild the historical membership from original Maxi pool FASTQs, follow the [historical reference guide](../../../docs/historical-maxi.md). Select `--maxi-method historical --maxi-reference DIR` for paired or panel analyses using that reference. Record the reference digest and the read mode with the results.
