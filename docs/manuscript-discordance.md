# Audit discordant manuscript calls

Use this audit after completing the 26-library manuscript comparison. It identifies which targets changed hit status, which threshold changed, and how the input counts compare. A hit requires both a donor-adjusted Δindel rate above 0.1 percentage points and a one-sided Benjamini–Hochberg q below 0.05. Each editor and strategy has its own 79-target correction.

## Run the audit

Start in the repository root with the analysis environment active. Set the paths to your completed run and the directory containing the same four submitted workbooks used for that run. The upstream checkout from the quick start supplies the public NoEP and HiFi allele tables. Python needs NumPy, SciPy and openpyxl.

```sh
python panel_report/audit_discordance.py \
  --run runs/manuscript_JOBID \
  --workbooks /path/to/submitted/workbooks \
  --upstream data/upstream \
  --output results/manuscript-discordance
```

Choose a new output directory. The command checks the workbook hashes against the run provenance and verifies the upstream Git commit. It uses the existing fitted models and library counts. If all 632 calls agree, it reports full agreement and exits. [Manuscript comparison workflow](../README.md#manuscript-comparison)

| Local output | What to inspect |
| --- | --- |
| `discordant_calls.tsv` | One row per discordant contrast, with the side that called a hit, the threshold crossed, fitted rates, Δ, p, q and observed or inferred counts for each arm. |
| `donor_counts.tsv` | Agent counts, public counts where available, and manuscript-compatible counts by donor and arm. The source column identifies how each event count was obtained. |
| `library_counts_check.tsv` | Counts by donor, replicate and arm, including the agent counts before cross-target filtering and matching public allele counts. |
| `reconstruction_qc.tsv` | Integer recovery and forward-fit errors, and whether the public inputs match the inferred manuscript margins. |

## Read a discrepancy

Start with `hit_in` and `gate`. The values `delta`, `q` and `delta+q` identify which threshold changed pass/fail status. Compare the adjusted Δ values in percentage points, then the two q values. The two additional q columns hold the focal p fixed while substituting the other 78 p values from the opposite source. They show whether the other targets alone could change the q decision.

Compare the indel numerator and total denominator separately for NoEP and the edited arm. Unfiltered counts reads. Dedup counts one representative allele per pBC family. The reported Δ comes from the donor-adjusted model. Pooling counts across donors gives a different estimand.

## Validation against the submitted tables

The Luminary refit agreed with the submitted tables on **608 of 632 calls**. Of the 24 discrepancies, 12 crossed only the Δ threshold, 11 crossed both thresholds and one crossed only the q threshold. Holding each focal p fixed and exchanging the other 78 p values preserved its q decision in all 24 cases.

| Editor and strategy | Only Δ crossed | Both crossed | Only q crossed |
| --- | ---: | ---: | ---: |
| WT Unfiltered | 6 | 3 | 0 |
| WT Dedup | 3 | 6 | 0 |
| WT Dedud+Dedup | 0 | 1 | 1 |
| HiFi Unfiltered | 2 | 0 | 0 |
| HiFi Dedup | 1 | 1 | 0 |

In all 18 WT Unfiltered and Dedup discrepancies, both the WT indel count and denominator differed from the counts inferred from the submitted model. Among the nine discordant WT Dedup cases, the agent counted 8.0%–14.1% as many WT families as the manuscript model implied. The corresponding NoEP denominators were 85.1%–92.7%. These ranges describe the discordant targets only.

## How to interpret inferred counts

The submitted workbooks report fitted probabilities and link-scale standard errors for nine donor-by-arm combinations. The script uses these to recover the binomial denominators, then uses the Firth adjusted-score equations to recover event totals by arm and by donor. It checks that the inferred counts round to integers and refits the coefficients from those integer margins. In the validation run it recovered the agent’s observed denominators and margins exactly for all 24 discordant models.

For the 22 Unfiltered and Dedup cases, direct public NoEP and HiFi counts agreed with the inferred donor denominators and arm event totals. The WT event count for each donor can then be estimated by subtracting those public events from the inferred donor total. For the two Dedud+Dedup cases, the public denominator margins differed, so the donor-by-arm event estimates are left blank.

The pinned public allele tables cover NoEP and HiFi. Assigning a manuscript WT difference to a specific replicate or read-processing step requires the original WT intermediate records. The public study separately provides WT FASTQs. [Upstream allele tables](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4/test/absolveseq_edits) · [Public FASTQ quick start](d3734-t1-quick-start.md)
