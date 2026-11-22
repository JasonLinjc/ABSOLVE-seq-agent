---
name: absolve-seq
description: Analyze original ABSOLVE-seq pooled reporter FASTQs end to end and deliver a report locally or on Slurm. Select a matched donor and replicate, compare plasmid barcode filtering criteria, export four indel-rate strategies, simulate detection power, deliver the report and figure, or reproduce manuscript treatment margins across the 79-target panel.
---

# ABSOLVE-seq analysis

Work from the repository root in the environment described in the README. Use the user's configured manifest and reference inputs. When neither a Slurm allocation nor `sbatch` is available, run the Python entry point directly on the local machine. Honor a request for local execution with `--executor local`.

For a public test dataset, follow the [D3734/T1 quick start](../../../docs/d3734-t1-quick-start.md#public-test-data-from-sra). Download both mates of NoEP `SRR39887646` and HiFi `SRR39887695` from study `SRP722226`. Use the example's `--ena` manifest option to verify the four files before submitting the pair with `--fresh`.

### Choices in chat

When a user asks to choose settings, offer short numbered options for the choices still open. Accept an option number or a plain-language answer. Wait for the answer before starting the selected analysis or rebuilding a reference. Carry choices forward across turns. For an explicit run with defaults, proceed once the inputs and sample identity are resolved.

- **Samples.** If the input source is open, offer the public D3734/T1 example or the user's original-study FASTQs and sample sheet. For user data, list the complete NoEP/HiFi pairs from the manifest. If several match, offer the actual donor/replicate combinations.
- **pBC annotation.** Offer three actions in this order: (1) use the published pBC annotation (default), (2) rebuild the historical annotation from the original study’s six Maxi FASTQs, or (3) reuse the historical annotation saved from a previous rebuild. For reuse, use the saved directory from the conversation or ask for its path. If a user offers a custom pool, explain that the builder accepts only the six original files listed in the [study file manifest](../../../runtime/maxi/inputs.sha256). This choice affects Dedud and Dedud+Dedup. Keep the selected read mode independent.
- **Execution.** Offer local execution or Slurm when both are available. For a remote cluster, resolve the SSH host and a shared project directory. Once the requested choices are settled, apply them to the commands below and continue through the report.

### Run and verify

1. For a donor/replicate request, run `python paired.py --list` and select the matching complete NoEP/HiFi pair. Run `python paired.py --donor DONOR --replicate REPLICATE` locally, or submit `sbatch paired.sbatch --donor DONOR --replicate REPLICATE` on Slurm. Use the requested Maxi method, with `manuscript` as the default.
2. For a manuscript-wide comparison, confirm that the four workbook extracts and 26 libraries are available. Follow the [paired-end reproduction guide](../../../docs/manuscript-reads.md) and run `python paired.py --manuscript-panel --read-mode manuscript` on the local machine or from a Slurm submission host. The run records the historical trimming steps for each library. To track changes against an earlier panel, supply its `panel_comparison.tsv` with `--baseline-panel`.
3. Check the local process exit code or Slurm job state, then verify the exported tables. For paired analyses verify 316 unique target/strategy rows in `results/paired_target_rates.tsv`. For the full panel verify 632 unique target/strategy/editor comparisons in the reported results directory’s `panel_comparison.tsv`.
4. Report rates in percent and differences in percentage points. For a pair, subtract NoEP from HiFi. For the manuscript panel, average the three donor-specific logits equally for each treatment before converting to rates. Report BH correction separately per strategy/editor. Distinguish threshold agreement from numerical agreement.
5. A fast matched-pair analysis may select `merged`, `any-pool`, `two-pools`, or `all-pools`. Those methods require the corresponding Maxi QC references. Use a distinct output directory for each run.

For the published annotation use `--maxi-method manuscript`. For a new build from the six original Maxi plasmid pool FASTQs, follow the [build instructions](../../../docs/historical-maxi.md#build-the-reference-once). For a newly rebuilt or reused historical annotation, pass its saved path as `--maxi-method historical --maxi-reference DIR`. Record the reference digest and the read mode with the results.


For a FASTQ-to-report request, carry the run through setup, sample matching, analysis, verification and delivery. Use the paths supplied by the user or the linked public example. Resolve the environment, references and manifest first. If donor or arm identity is ambiguous, ask the user which comparison to run. Once the inputs are resolved, continue through the remaining steps in the same conversation.

For Slurm, wait for the final aggregation job to succeed. Run `python verify_results.py power RESULTS` and check the linked outputs in `report.md`, including `power_estimates.tsv` and `power_provenance.json`. Describe the assumed effects, simulation count and nominal significance level. Report a missing denominator as unavailable. Return the report, effect and power tables, and links to the Dedud PDF and Excel. The pair produces a single-donor figure. The panel produces Figure 3b with newly simulated power. For a matching study benchmark, run the appropriate `verify_results.py` profile. If the user explicitly requests a historical redraw using submitted power, follow the [figure skill](../absolve-seq-figures/SKILL.md). On failure, identify the stage and log.
