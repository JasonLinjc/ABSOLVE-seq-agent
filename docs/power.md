# Estimating editing effects and detection power

Run the ordinary paired workflow or full manuscript panel. Both generate effects from the input counts, simulate power and write a report, figure and source data into the same results directory. The count refit uses the same final steps.

| Run | Editing effect | Power model |
| --- | --- | --- |
| Matched NoEP and HiFi libraries | Difference between observed indel fractions, in percentage points | Two-group binomial design at the observed arm depths. The control baseline is the Firth estimate `(NoEP indels + 0.5) / (NoEP denominator + 1)`. |
| Full panel | Model-based WT and HiFi rates minus NoEP | Preserve the observed donor-by-arm depths and each donor's fitted NoEP probability. Simulate all three arms jointly and refit an additive treatment-and-donor Firth logistic regression for each draw. |

For each target and strategy, the agent adds 0.01 or 0.1 **percentage points** to the control probability in each edited arm. It simulates both control and treated counts using binomial sampling. Within a donor and arm, summing replicate depths is equivalent to summing independent binomial draws with the same underlying probability. The manuscript mode uses equal donor weights on the logit scale for its marginal editing estimates. The power test uses the fitted editor coefficient and its model-based standard error. [Original manuscript power script](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/absolveseq/8_power_analysis.R)

**Interpretation.** Power is the proportion of 1,000 simulations that detect a positive editor effect at one-sided nominal *p* < 0.05. The observed panel calls additionally use BH-adjusted q-values and an effect-size cutoff. Simulations condition on the measured depths, fitted background and binomial independence of the counted units. For a single pair they describe sampling variability within that pair. Generalization across biological donors requires a replicated design. A zero denominator or infeasible assumed effect has a blank estimate and an explicit status.

**Files.** Open `report.md` for the editing and power summaries. Each row of `power_estimates.tsv` identifies the target, strategy, editor and assumed effect, along with detections, simulation count, power and a 95% Wilson interval for Monte Carlo error. The fixed base seed and input hashes are in `power_provenance.json`. Changing inputs can change the estimates. Each target, strategy and effect has a separate deterministic random stream.

The PDF, PNG and Excel table under `figures/` show Dedud. The paired plot orders targets by ID. The full panel uses the paper's Figure 3b target order and exports 158 target/editor rows. Historical redraws can import the original submitted power values explicitly. [Historical Figure 3b redraw](manuscript-figures.md)

## Verify a completed run

Run `python verify_results.py power PATH/TO/RESULTS`. The check verifies input hashes, coverage, simulation counts, interval bounds, and the plotted effects and power against the run tables and workbook. The existing effect-table benchmarks and historical Figure 3b redraw profiles remain available separately.

## Study validation

For the full panel, the 8,216 archived count rows reproduced the accepted effect-table checksum and all 632 submitted primary calls. Across the 632 observed editor coefficients, the largest absolute difference between the simulation solver's Wald z-statistics and the independent R/brglm2 fit was 0.0000034. The power TSV has 1,264 rows for the two effect sizes.

For the 158 Dedud target/editor combinations in Figure 3b, the median absolute differences between independently simulated and submitted power were 0.35 percentage points at an assumed editing effect of 0.01 pp and 0 percentage points at 0.1 pp. The corresponding maximum differences were 5.1 and 0.4 percentage points. These differences compare detection probabilities from separate simulations.
