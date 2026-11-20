---
name: absolve-seq-figures
description: Regenerate ABSOLVE-seq manuscript panels 2a, 2b, 2d, 3b and Supplementary Figures 10–13 from the original 79-target study. Use for the published editing estimates, confidence intervals, threshold calls, and candidate-specific power panels.
---

# Reproduce manuscript figures

Use the submitted manuscript figure identifier to choose the analysis. Work from the repository root with Python, Matplotlib, NumPy and openpyxl. The source workbooks are the original `SupplementaryTable2.xlsx` through `SupplementaryTable5.xlsx`.

| Figure | Strategy | Source |
|---|---|---|
| 2a | All four | Archived allele tables grouped by target and arm |
| 2b | All four | Historical NoEP_read_melt.csv |
| 2d | Maxi yield | Published per-target index built by setup_references.py |
| 3b | Dedud | Supplementary Table 2 |
| S10 | Unfiltered | Supplementary Table 3 |
| S11 | Dedup | Supplementary Table 4 |
| S12 | Dedud+Dedup | Supplementary Table 5 |
| S13a,b | Four strategies | Original draft5 plotting tables for the artwork or Supplementary Tables 2–5 for the final-workbook view |

Use `scripts/figure2.py` for the coverage, background and Maxi panels. Use `scripts/figure2_coverage.py ARCHIVE > coverage.tsv` to obtain the 948 input measurements for 2a. The archive layout and source-specific commands are in the [figure guide](../../../docs/manuscript-figures.md).

For a redraw of the submitted estimates, supply `--tables`. For a redraw from a completed agent model, also supply `--results` containing the matching `panel_comparison.tsv` and `model_field_comparison.tsv`. Both workbook and refit modes use the candidate-specific power in the submitted workbooks. S13 also accepts the original draft5 summary and power tables with `--figure-tables`. Use that mode for a reproduction of its printed bar heights. A two-library paired run cannot supply the three-donor model or the archived power simulations.

```bash
plot=.agents/skills/absolve-seq-figures/scripts/manuscript_figures.py
python "$plot" 3b --tables submitted-workbooks --results runs/panel-refit --out figures/3b
python "$plot" S13 --figure-tables draft5/tables --out figures/S13
```

Figures 3b and S10–12 split sgRNAs 1450 and 1617, following the [submitted target order](references/target_order.tsv). Each treatment confidence limit is shifted by the NoEP point estimate, as in the original plotting code. The five heatmap rows show effect thresholds, adjusted p-values and candidate-specific power. HiFi occupies the upper-left triangle and SpCas9 the lower-right. Figure S13 uses the combined effect and BH thresholds on the left and plots each candidate's power on the right following the submitted artwork.

Open the PDF and inspect both guide blocks and the heatmap. Verify 158 plotted target/editor rows per single-strategy figure and 632 for S13. Return the PDF, plotted data and provenance. State whether estimates are from the submitted tables or the agent refit and that power comes from the submitted simulations.

For S13, identify the input version explicitly. The WT Dedud count at Δ > 0.01 percentage points and BH q < 0.05 is seven in the artwork and eight from the final workbook.

For other panels, identify which original inputs they need before promising a redraw. Figure 2c needs indel-level pBC quality classes, Figure 3a,c require endogenous measurements as well as reporter estimates, and Figure 4/S14 require outcome-level Lindel predictions and experimental pBC support. A saved artwork or an illustrative diagnostic is a separate deliverable from a data-driven redraw.
