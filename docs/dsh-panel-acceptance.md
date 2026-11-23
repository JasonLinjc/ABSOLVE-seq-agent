# DeepSeek Harness acceptance from raw FASTQs to a full-panel report

On November 23, 2026, DeepSeek Harness (`dsh` 0.2.2-rc.1, model `DeepSeek-V4-Flash-0731`) submitted a fresh run on Luminary through a site-configured MCP tool. The run processed 52 staged FASTQs covering all 79 targets in 26 libraries from three donors (9 NoEP, 9 WT and 8 HiFi). It used uniform paired-end CRISPResso2 2.2.13 processing and the published plasmid barcode annotation. [Tested agent revision](https://github.com/JasonLinjc/ABSOLVE-seq-agent/tree/e15d0e73d5f8313742fa8ad8be3d54baf3af6481)

## Execution and verification

The submission tool invoked the following command in a new run directory. Slurm completed 26 preparation tasks, 26 alignment tasks and the final report job in 14 minutes from the first task start. The same `dsh` session queried the scheduler and received `PASS`. The statistical tables, report links and figure workbook were checked locally after download. [Full-panel inputs and setup](manuscript-reads.md#analyze-all-study-samples)

```bash
python paired.py --manuscript-panel --read-mode manuscript --fresh --output run
```

| Check | Result |
| --- | --- |
| Raw input and target coverage | The trimming logs recorded 334,970,765 input read pairs. All 2,054 library/target combinations had fresh count records. |
| Editing estimates | All 8,216 library count rows and 632 modeled contrasts were present. Both TSV files matched the previous uniform paired-end run byte for byte. |
| Power and Figure 3b | `python verify_results.py power run/results` passed for 1,264 power rows and all 158 plotted target/editor combinations, including their Excel values. |
| Download integrity | Local SHA-256 digests matched the Luminary digests for all 16 result files. |

## Editing calls

Under Dedud, all three positive controls passed for both WT and HiFi versus NoEP (Δ > 0.1 percentage points and BH-adjusted q < 0.05). These were the on-target controls `1450_OT_0000_REF` and `1617_OT_0000_REF` and the off-target control `1617_OT_0040_ALT`. WT also passed at `1450_OT_0071_REF` (Δ = 0.147 percentage points, q = 1.46 × 10⁻¹⁵).

Across all four strategies, 628 of the 632 calls matched the submitted supplementary tables. The four differences occurred in WT (two Unfiltered and two Dedup). The separate historical-processing run reproduced all 632 submitted calls. [Historical panel acceptance](historical-panel-acceptance.md)

## SHA-256 for the two archive comparisons

| File | Digest |
| --- | --- |
| `library_counts.tsv` | `aefad456069a4eb3d8eef4e2cdef9b913c6e8d2791df2be04ea306f60b9fef40` |
| `panel_comparison.tsv` | `3f52e3068390b9bf5294a816180d2dde7107a9cd29458f0f50b0079341c6290d` |
