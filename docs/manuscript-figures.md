# Regenerate the submitted ABSOLVE-seq figures

The figure identifiers refer to the original 79-target manuscript submitted in September 2026. The [figure skill](../.agents/skills/absolve-seq-figures/SKILL.md) writes a vector PDF, a PNG, the plotted data and the SHA-256 of the inputs. Use a new output directory for each run. Execute these examples from the repository root.

## Coverage, background and Maxi yield

```bash
scripts=.agents/skills/absolve-seq-figures/scripts
python "$scripts/figure2_coverage.py" /archive/ABSOLVE-seq > coverage.tsv
python "$scripts/figure2.py" 2a coverage.tsv --out figures/2a
python "$scripts/figure2.py" 2b /archive/ABSOLVE-seq/NoEP_read_melt.csv --out figures/2b
python "$scripts/figure2.py" 2d data/references/published --out figures/2d
```

- **Figure 2a.** The archive contains `raw` and `dedudS_filtered_v3` directories of per-replicate `*annot.txt` files under `target_umi_barcode_table_WTCas9_NoEP_filtered` and `target_umi_barcode_table_HiFi_NoEP`. Pool within each target and treatment. Sum supporting reads for Unfiltered and Dedud. Count distinct pBC identifiers across the treatment for Dedup and Dedud+Dedup. Summing per-library family counts can count a recaptured barcode multiple times. The 79 targets, three treatments and four strategies yield 948 plotted values.
- **Figure 2b.** Read the 316 target/strategy percentages from the archived NoEP plotting CSV. Use the original `value` column on the linear axis. The unfiltered range is 0.0532%–0.8480%. The arithmetic mean of the archived coordinates is 0.1884%. The manuscript text reports 0.16%.
- **Figure 2d.** Run the standard `setup_references.py` on the pinned upstream checkout to generate the published index. The figure uses 11,239,298 unfiltered and 5,190,216 high-quality target/pBC records. Both counts were independently checked per target against the original Maxi archive. The historical analysis reference additionally retains 58,725 low-entropy two-tBC memberships. Using that reference for the figure would change the right-hand distribution.

The small displacement of overlapping points is for readability. The two on-target controls and the CPS1 positive control have red triangle markers. Figure 2c remains outside the verified plotting commands. A replay from archived overlapping pBC quality classes yielded different percentages from the submitted coordinates.

## Editing and candidate-specific power

```bash
plot=.agents/skills/absolve-seq-figures/scripts/manuscript_figures.py
python "$plot" 3b --tables submitted-workbooks --out figures/3b-submitted
python "$plot" 3b --tables submitted-workbooks --results runs/panel-refit --out figures/3b-refit
python "$plot" S10 --tables submitted-workbooks --out figures/S10
python "$plot" S11 --tables submitted-workbooks --out figures/S11
python "$plot" S12 --tables submitted-workbooks --out figures/S12
```

The four figures use Dedud, Unfiltered, Dedup and Dedud+Dedup respectively. Each has 79 targets and two editor estimates. The target order follows the submitted Figure 3b artwork, with label identities checked against Supplementary Table 1. Confidence bounds are the editor's marginal bounds shifted by the unedited estimate, following the historical plotting code. Heatmaps combine Δ thresholds, within-editor/strategy BH adjustment and the corresponding workbook's simulated power. With `--results`, effect estimates, interval bounds and q-values come from the agent's equal-donor model output. The workbook supplies the power columns. No simulation is repeated during plotting.

## Supplementary Figure 13 and its input versions

```bash
python "$plot" S13 --figure-tables /archive/tablesAndFigures_draft5/tables --out figures/S13-original
python "$plot" S13 --tables submitted-workbooks --out figures/S13-workbooks
python "$plot" S13 --tables submitted-workbooks --results runs/panel-refit --out figures/S13-refit
```

To reproduce the original artwork, supply its four `summaryResult_*_fdr.tsv` files and eight `powerAnalysis_*_0.01.csv` or `*_0.1.csv` files. The original plotting notebook was `tablesAndFigures_draft5/visualize.ipynb`. The left panels count targets above the effect threshold with BH q below 0.05. The right panels show individual target power values and reference lines at 0.8 and 0.9, as in the submitted artwork.

The original figure's WT Dedud bar has seven targets at Δ above 0.01 percentage points. The final workbook and its refit give eight. For `1450_OT_0210_ALT`, the original plot table gives Δ = 0.00847315 percentage points and the submitted workbook gives 0.01050316. Both q-values are below 0.05. At the manuscript's primary 0.1 percentage-point threshold, all 632 panel classifications agree.

Figures 1 and S9 are schematics. Figure 3a,c includes endogenous measurements and Figure 4/S14 includes outcome-level predictions and family support. Those panels have different source requirements. The commands above regenerate the listed quantitative panels.
