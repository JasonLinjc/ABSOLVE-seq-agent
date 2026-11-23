---
name: absolve-seq
description: Analyze original ABSOLVE-seq pooled reporter FASTQs end to end and deliver a report locally or on Slurm. Select a matched donor and replicate, compare plasmid barcode filtering criteria, export four indel-rate strategies, simulate detection power, deliver the report and figure, or reproduce manuscript treatment margins across the 79-target panel. Use this skill for requests in any language to validate pBC annotation reproducibility on Luminary or audit historical pBC family completeness across WT, NoEP and HiFi.
---

# ABSOLVE-seq analysis

Work from the repository root in the environment described in the README. Use the user's configured manifest and reference inputs. When neither a Slurm allocation nor `sbatch` is available, run the Python entry point directly on the local machine. Honor a request for local execution with `--executor local`.

For a public test dataset, follow the [D3734/T1 quick start](../../../docs/d3734-t1-quick-start.md#public-test-data-from-sra). Download both mates of NoEP `SRR39887646` and HiFi `SRR39887695` from study `SRP722226`. Use the example's `--ena` manifest option to verify the four files before submitting the pair with `--fresh`.

For a historical per-pBC CRISPRessoBatch directory, use the [historical allele collector](../../../docs/manuscript-reads.md#collect-a-historical-per-pbc-batch). Start downstream counting only after the collector succeeds. If collection fails, report the named pBCs and affected tables.

### Validate pBC annotation reproducibility

Route “Validate pBC annotation reproducibility” and equivalent requests in any language directly here, including a short Chinese request with no host, donor or mode specified. When offering tasks, list this as an optional full-panel check. The request selects Luminary, all 26 libraries and both read modes. It compares the published annotation with the saved six-Maxi rebuild under the published rule using existing alignments.

1. From the client checkout run `python3 luminary_pbc.py`. The launcher uses the configured Luminary checkout and prints the job ID and local delivery directory. If the user supplies a checkout or SSH alias, pass `--repo` or `--host`. For a missing configuration, follow the [one-time setup](../../../docs/pbc-validation.md#run-from-codex-or-claude-code). Skip donor selection, annotation rebuilding and local execution for this request.
2. Follow that process through completion. If monitoring or delivery was interrupted, run `python3 luminary_pbc.py --resume LOCAL_DELIVERY_DIR` with the directory from the same request. If submission was interrupted before a job ID was recorded, check the Luminary queue and accounting before attempting another submission.
3. Read `summary.tsv` and `differences.tsv` from the verified delivery. In the user's language, give the job ID, counts and model rows checked per mode and observed differences. Link the local `summary.tsv`, `differences.tsv` and `SHA256SUMS`. Exit 1 delivers completed comparisons with differences. For exit 2, identify the failed stage using the job receipt and Slurm log when available.

### Audit historical pBC family completeness

For a request to check missing families across WT, NoEP and HiFi, run the [full-panel family acceptance command](../../../docs/historical-panel-acceptance.md#rerun-the-full-panel-pbc-family-audit) with the historical project archive and the receipt-checked FASTQ family replays. Use `--expect-study` to verify the documented study result. Return the batch summary and per-pBC differences with their evidence source and the changes in both filtered numerators and denominators. A matching historical baseline includes the documented omission. A request to compare Maxi annotation methods uses the separate annotation validation command above.

### Choices in chat

When a user asks to choose settings, offer short numbered options for the choices still open. Accept an option number or a plain-language answer. Wait for the answer before starting the selected analysis or rebuilding a reference. Carry choices forward across turns. For an explicit run with defaults, proceed once the inputs and sample identity are resolved.

- **Samples.** If the input source is open, offer the public D3734/T1 example or the user's original-study FASTQs and sample sheet. For user data, list the complete NoEP/HiFi pairs from the manifest. If several match, offer the actual donor/replicate combinations.
- **pBC annotation.** Offer three actions in this order: (1) use the published pBC annotation (default), (2) rebuild from the original study’s six Maxi FASTQs using the published entropy rule, with the historical two-tBC exception available on request, or (3) reuse a saved rebuild. For reuse, use the saved directory from the conversation or ask for its path. Read its `provenance.json` and select the recorded method. If the user explicitly asks for a different rule, explain the mismatch and offer a matching rebuild. If a user offers a custom pool, explain that the builder accepts only the six original files listed in the [study file manifest](../../../runtime/maxi/inputs.sha256). This choice affects Dedud and Dedud+Dedup. Keep the selected read mode independent.
- **Execution.** Offer local execution or Slurm when both are available. For a remote cluster, resolve the SSH host and a shared project directory. Once the requested choices are settled, apply them to the commands below and continue through the report.

### Run and verify

1. For a donor/replicate request, run `python paired.py --list` and select the matching complete NoEP/HiFi pair. Run `python paired.py --donor DONOR --replicate REPLICATE` locally, or submit `sbatch paired.sbatch --donor DONOR --replicate REPLICATE` on Slurm. Use the requested Maxi method, with `manuscript` as the default.
2. For reproduction of the submitted results from original FASTQs, follow the [historical full-panel guide](../../../docs/manuscript-reads.md#reproduce-the-submitted-calls-from-the-original-fastqs). Select `--manuscript-panel --read-mode historical --fresh`, process all 26 libraries and wait for the finish job. Run the `panel-historical` workbook check and the `power` check, then link the report, comparison table, Figure 3b and Excel. For the uniform CRISPResso2 2.2.13 manuscript workflow, confirm that the four workbook extracts and 26 libraries are available. Follow the [paired-end reproduction guide](../../../docs/manuscript-reads.md) and run `python paired.py --manuscript-panel --read-mode manuscript` on the local machine or from a Slurm submission host. The run records the historical trimming steps for each library. To track changes against an earlier panel, supply its `panel_comparison.tsv` with `--baseline-panel`.
3. Check the local process exit code or Slurm job state, then verify the exported tables. For paired analyses verify 316 unique target/strategy rows in `results/paired_target_rates.tsv`. For the full panel verify 632 unique target/strategy/editor comparisons in the reported results directory’s `panel_comparison.tsv`.
4. Report rates in percent and differences in percentage points. For a pair, subtract NoEP from HiFi. For the manuscript panel, average the three donor-specific logits equally for each treatment before converting to rates. Report BH correction separately per strategy/editor. Distinguish threshold agreement from numerical agreement.
5. A fast matched-pair analysis may select `merged`, `any-pool`, `two-pools`, or `all-pools`. Those methods require the corresponding Maxi QC references. Use a distinct output directory for each run.

For the published annotation use `--maxi-method manuscript`. For a new build from the six original Maxi plasmid pool FASTQs, follow the [build instructions](../../../docs/historical-maxi.md#build-the-reference-once). Use `--published-maxi` to rebuild the published rule or `--historical-maxi` for the historical rule. Pass the saved path to the analysis with `--maxi-reference DIR` and the corresponding `--maxi-method manuscript` or `--maxi-method historical`. Record the reference digest and the read mode with the results.


For a FASTQ-to-report request, carry the run through setup, sample matching, analysis, verification and delivery. Use the paths supplied by the user or the linked public example. Resolve the environment, references and manifest first. If donor or arm identity is ambiguous, ask the user which comparison to run. Once the inputs are resolved, continue through the remaining steps in the same conversation.

For Slurm, wait for the final aggregation job to succeed. Run `python verify_results.py power RESULTS` and check the linked outputs in `report.md`, including `power_estimates.tsv` and `power_provenance.json`. Describe the assumed effects, simulation count and nominal significance level. Report a missing denominator as unavailable. Return the report, effect and power tables, and links to the Dedud PDF and Excel. The pair produces a single-donor figure. The panel produces Figure 3b with newly simulated power. For a matching study benchmark, run the appropriate `verify_results.py` profile. If the user explicitly requests a historical redraw using submitted power, follow the [figure skill](../absolve-seq-figures/SKILL.md). On failure, identify the stage and log.
