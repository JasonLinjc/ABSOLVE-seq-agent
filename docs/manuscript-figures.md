# Reproduce Figure 3b and its source data

Figure 3b shows the Dedud analysis for 79 targets and two editors, WT SpCas9 and HiFi, each compared with NoEP. Run from the repository root with the [locked analysis environment](container.md). Both commands require the original submitted `SupplementaryTable2.xlsx` in the directory passed to `--tables`.

On Linux x86-64, install the exact environment from this checkout and run its scripts. If using an older released container, bind this checkout and run the plotting script from the mounted repository path.

```bash
conda create --yes --name absolve-figure3b --file docs/environment-linux-64.lock
conda activate absolve-figure3b
```

## Draw the submitted estimates

```bash
plot=.agents/skills/absolve-seq-figures/scripts/manuscript_figures.py
python "$plot" 3b --tables submitted-workbooks --out figures/3b-submitted
python verify_results.py figure3b-submitted figures/3b-submitted
```

## Draw the agent refit

Refit the [26-library count table](historical-maxi.md#refit-existing-counts) or supply the directory containing its `panel_comparison.tsv` and `model_field_comparison.tsv`. Use a new output directory for each draw.

```bash
plot=.agents/skills/absolve-seq-figures/scripts/manuscript_figures.py
python "$plot" 3b --tables submitted-workbooks --results runs/panel-refit --out figures/3b-refit
python verify_results.py figure3b-refit figures/3b-refit
```

Each command writes `Figure3b.pdf`, `Figure3b.png`, `Figure3b_source_data.xlsx`, `plotted_data.tsv` and `provenance.json`. The first worksheet contains the exact plot records at Excel numeric precision. The Notes sheet identifies units, thresholds, estimate source and input hashes. Each row pairs one target with one editor. The target order repeats for each editor.

Δ is the edited marginal rate minus the NoEP marginal rate in percentage points. Each displayed interval shifts the treatment confidence bounds by the NoEP point estimate. The heatmap applies strict Δ thresholds of 0.01 and 0.1 percentage points and a within-editor BH q threshold of 0.05. Both modes use target-specific simulation power from Supplementary Table 2. HiFi occupies the upper-left triangle and WT the lower-right.

The benchmark verifies the input fingerprints, a SHA-256 of all 158 ordered TSV records, and every Excel data cell. Numeric Excel values allow a relative tolerance of 2 × 10⁻¹⁵ for serialization, with exact zeros. The two estimate sources have separate reference digests. Image rendering can vary with fonts and graphics libraries.
