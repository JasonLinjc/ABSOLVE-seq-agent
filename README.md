<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/absolve-seq-agent-dark.png" />
    <img src="docs/assets/absolve-seq-agent.png" width="480" alt="ABSOLVE-seq Agent" />
  </picture>
</p>

<h1 align="center">Quantify off-target gene editing with ABSOLVE-seq</h1>

<p align="center">
  A conversational analysis agent for the original 79-target reporter panel, available in Codex and Claude Code.
</p>

<p align="center">
  <a href="#start-with-one-message">💬 <strong>Start in chat</strong></a> &nbsp;·&nbsp;
  <a href="docs/d3734-t1-quick-start.md">🧬 Public dataset tutorial</a> &nbsp;·&nbsp;
  <a href="#study-panel-results">📊 Read the results</a> &nbsp;·&nbsp;
  <a href="#figure-3b-and-reproducibility">📄 Manuscript Figure 3b</a>
</p>

**What is ABSOLVE-seq?** ABSOLVE-seq (Assessment By Stand-in Off-target LentiViral Ensemble with sequencing) measures insertions and deletions (indels) in reporters carrying reference or variant sequences from candidate CRISPR off-target sites. The agent estimates how much editing differs between each Cas9 treatment (WT or HiFi) and the NoEP control at each target, and how likely the experiment is to detect a small increase. [Original study](https://doi.org/10.1101/2024.07.24.605019)

## Example conversations

<a id="from-paired-fastqs-to-a-report"></a>

### Start with one message

Open a folder in Codex or Claude Code and send the **User** message below. The agent can use your study manifest or download the 52 public FASTQs. The current panel workflow also needs the model references extracted from Supplementary Tables 2–5, or the four workbooks to prepare them. [Full-panel setup](docs/manuscript-reads.md#analyze-all-study-samples)

<table>
<tbody>
<tr><td><strong>👤 User · One-message prompt</strong><br>Set up https://github.com/JasonLinjc/ABSOLVE-seq-agent and analyze all 79 targets in all 26 original study libraries from raw FASTQs. Use my configured data or fetch the public reads. Include every donor and replicate. For each of the four filtering strategies, compare WT and HiFi separately with NoEP. Report indel rates, differences in percentage points, p-values, adjusted q-values and detection power. Verify the outputs, summarize the findings here, and link the report, complete tables, Figure 3b and its Excel source data.</td></tr>
<tr><td><strong>🤖 Agent · Example result</strong><br>Across 9 NoEP, 9 WT and 8 HiFi libraries from three donors, all three positive controls passed the Dedud thresholds (Δ &gt; 0.1 percentage points and q &lt; 0.05) for both editors versus NoEP. These were the two on-target controls, <code>1450_OT_0000_REF</code> and <code>1617_OT_0000_REF</code>, and the previously validated off-target control, <code>1617_OT_0040_ALT</code>.<br><br>At the off-target control, the estimated indel rate was 0.01% for NoEP, 51.77% for WT and 4.49% for HiFi. One additional target, <code>1450_OT_0071_REF</code>, passed the Dedud thresholds for WT only.<br><br>Open <code>runs/study-panel/results/report.md</code> for the results across all targets and strategies, including detection power and the figure source data.</td></tr>
</tbody>
</table>

<details>
<summary><strong>Try a smaller example or select one matched pair</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Run the public D3734/T1 HiFi and NoEP example from FASTQs and show me the report.</td></tr>
<tr><td><strong>👤 User · With local data</strong><br>My FASTQs are in <code>/path/to/fastqs</code> and my sample sheet is at <code>/path/to/sample-sheet</code>. Find the matched HiFi and NoEP pairs. If several are available, let me choose one, then run it through to the report.</td></tr>
</tbody>
</table>

The public pair uses four FASTQs (1.58 GB). [Paired tutorial](docs/d3734-t1-quick-start.md) · [Docker acceptance record](docs/container.md#d3734t1-docker-acceptance)

For remote execution, include the cluster’s SSH alias and a shared working directory. Give paths that the compute nodes can read. [Cluster setup and transfer](docs/container.md#stage-public-inputs-from-a-networked-machine)

</details>

<a id="choose-settings-in-chat"></a>
<details>
<summary><strong>Choose a plasmid barcode annotation</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>I generated a pBC annotation from the study’s plasmid pools last week. Can I reuse it for this pair?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Yes. Give me the saved annotation directory and I’ll use the rule recorded there. You can also use the published annotation or ask me to rebuild it from the original study’s six Maxi FASTQs using the published entropy threshold or the historical two-tBC rule.</td></tr>
</tbody>
</table>

[Maxi inputs and build instructions](docs/historical-maxi.md#build-the-reference-once)

</details>

<a id="manuscript-comparison"></a>
<details>
<summary><strong>Reproduce Supplementary Tables 2–5</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>I have the original 26 libraries and the submitted workbooks. Can you reproduce Supplementary Tables 2–5?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>I’ll use the historical processing for WT, NoEP and HiFi and the historical pBC annotation. I’ll compare the indel rates, editing differences, q-values and positive calls with the workbooks across all 632 contrasts.</td></tr>
</tbody>
</table>

[Historical runtime and inputs](docs/manuscript-reads.md#reproduce-the-submitted-calls-from-the-original-fastqs) · [Full-panel acceptance](docs/historical-panel-acceptance.md)

</details>

<a id="validate-pbc-annotation-reproducibility"></a>
<details>
<summary><strong>Validate pBC annotation reproducibility</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Validate pBC annotation reproducibility.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>I’ll compare the published and rebuilt annotations using the saved alignments for all 26 libraries in the configured Luminary archive. I’ll return the summary, field-level differences and SHA-256 manifest.</td></tr>
</tbody>
</table>

[Archive setup and output guide](docs/pbc-validation.md#run-from-codex-or-claude-code)

</details>

## How the agent works

<a id="fastq-to-report-workflow"></a>

### Agent architecture

<p align="center">
  <img src="docs/assets/agent-architecture.svg" width="660" alt="Codex and Claude Code use shared skills to run the analysis locally or on Slurm. The workflow connects sequencing and barcode inputs to editing and power estimates, a report, a figure and Excel data." />
</p>

**Inputs.** Each library contributes an R1 and an R2 FASTQ. The manifest maps both files to the sample and records their SHA-256 checksums. [Manifest template](examples/manifest.example.tsv)

**Processing.** Trimmomatic 0.36 trims both mates. The target barcode (tBC) assigns each read pair to a target. CRISPResso2 2.2.13 uses FLASH to merge the mates, aligns the merged reads to the target amplicon and calls indels in the cut-site window. Reads from the same sample and target are grouped by plasmid barcode (pBC) for counting.

## Study panel results

Open `results/report.md` in the run directory. The model accounts for donor and gives the three donors equal weight when estimating each treatment’s indel rate.

| File in `results/` | What it contains |
| --- | --- |
| `library_counts.tsv` | 8,216 indel counts and denominators (26 libraries × 79 targets × 4 strategies). |
| `panel_comparison.tsv` | 632 contrasts (WT or HiFi versus NoEP). Each row gives the fitted rates, Δ, one-sided p-value, and Benjamini–Hochberg q-value adjusted across 79 targets within that editor and strategy. |
| `power_estimates.tsv` | 1,264 estimates of detection power, one for each contrast at assumed increases of 0.01 and 0.1 percentage points. |
| `figures/Figure3b.pdf` and `figures/Figure3b_source_data.xlsx` | Dedud results for both editors, with 158 target/editor combinations in the figure and workbook. |

## Paired results

Open the `results/` folder inside your run directory.

| File | What it tells you |
| --- | --- |
| `report.md` | Start here for the summary and links to each output. |
| `paired_target_rates.tsv` | Indel counts, denominators and rates (%) for both arms. Δindel rate is HiFi minus NoEP in **percentage points**. A zero denominator leaves the rate undefined. |
| `power_estimates.tsv` | Estimated probability of detecting an increase of 0.01 or 0.1 percentage points at the observed depth for each target and strategy. [Power method](docs/power.md) |
| `figures/paired_effect_power.{pdf,png,xlsx}` | The Dedud figure for this pair and the Excel table containing the values plotted. |

The four strategies below use the same alignments. By default, the published plasmid barcode annotation determines which pBCs qualify for each target.

| Strategy | What counts toward the indel rate |
| --- | --- |
| Unfiltered | All aligned reads |
| Dedup | One representative allele per sample, target and pBC |
| Dedud | Reads with a qualifying pBC that pass cross-target filtering |
| Dedud+Dedup | One representative allele per retained pBC family |

The paired workflow selects the most frequent allele in each family for deduplication. It breaks ties by aligned read sequence, then by aligned reference sequence. Library counts, sample identity and processing settings are recorded in `library_counts.tsv`, `sample_inventory.tsv` and `provenance.json`.

## Set up in Codex or Claude Code

For manual setup, clone the repository and open the checkout in your signed-in client.

```bash
git clone https://github.com/JasonLinjc/ABSOLVE-seq-agent.git
cd ABSOLVE-seq-agent
```

| Client | How to start |
| --- | --- |
| Codex | Open the folder in the Codex app or run `codex` in the checkout. Ask for the ABSOLVE-seq skill or type `$absolve-seq` in the CLI. [Codex quick start](https://developers.openai.com/codex/quickstart/) |
| Claude Code | Run `claude` in the checkout and enter `/absolve-seq` followed by your request. [Claude Code quick start](https://code.claude.com/docs/en/quickstart) |

For the analysis runtime, download the prebuilt Linux x86-64 Docker or Apptainer image, or install the locked Conda environment. [Download v0.2.0](https://github.com/JasonLinjc/ABSOLVE-seq-agent/releases/tag/v0.2.0) · [Load and run the image](docs/container.md#download-a-prebuilt-image) · [Conda setup](docs/d3734-t1-quick-start.md#1-get-the-agent-and-create-an-environment)

<details>
<summary><strong>Check sample discovery before downloading data</strong></summary>

```bash
ABSOLVE_MANIFEST=examples/manifest.example.tsv python paired.py --list
```

The example manifest lists `NoEP_pool_D3734_T1` and `HiFi_pool_D3734_T1`. Its FASTQ paths are placeholders.

</details>

## Quick start

For all 26 libraries, follow the [full-panel walkthrough](docs/manuscript-reads.md#analyze-all-study-samples).

For a terminal run, follow the [D3734/T1 walkthrough](docs/d3734-t1-quick-start.md) to prepare the reference index and manifest and verify the NoEP (`SRR39887646`) and HiFi (`SRR39887695`) FASTQs. Run the following in the paired environment after preparing those inputs.

```bash
python paired.py --executor local --donor D3734 --replicate T1 --threads 4 --output runs/D3734-T1
python verify_results.py pair-manuscript runs/D3734-T1/results
python verify_results.py power runs/D3734-T1/results
```

The second command checks the public pair against its accepted count and rate tables. The third checks the simulated power estimates, figure data and Excel values. Use a new output directory for each analysis.

On a Slurm submission host, omit `--executor local` and `--threads 4`. The agent submits preparation, alignment and reporting jobs. Wait for the finish job recorded in `jobs.json` before opening the report. [Slurm instructions](docs/d3734-t1-quick-start.md#slurm-submission)

<a id="local-inputs-and-caches"></a>
<details>
<summary><strong>Configure local inputs</strong></summary>

`ABSOLVE_DATA` selects the data directory (default `data/` in the checkout). `ABSOLVE_MANIFEST` defaults to `$ABSOLVE_DATA/sample_manifest.tsv`. For the original laboratory files, `python -m tools.resolve_inputs INVENTORY.tsv METADATA.xlsx OUTPUT.tsv` joins both mate filenames from the checksum-verified inventory to the submission workbook’s `SRA_data` sheet.

</details>

## Figure 3b and reproducibility

The full panel combines three donors to generate `Figure3b.pdf` and `Figure3b_source_data.xlsx`, with 158 plotted target/editor rows. The default panel figure uses newly simulated power. To redraw the published graphic with its original power values, ask for **“Figure 3b with submitted power.”** [Figure 3b commands and verification](docs/manuscript-figures.md)

## Repository layout

The three main commands live at the top level. Use [setup_references.py](setup_references.py) to prepare the reference data, [paired.py](paired.py) to run the analysis, and [verify_results.py](verify_results.py) to check the results.

| Directory | What you’ll find |
| --- | --- |
| [examples/](examples/) and [docs/](docs/) | Sample manifests, public FASTQ links, walkthroughs and acceptance records |
| [.agents/skills/](.agents/skills/) and [.claude/skills/](.claude/skills/) | Analysis and figure instructions for Codex and Claude Code |
| [absolve/](absolve/) | Barcode filtering, target identity, power estimation, figures and reports |
| [panel_report/](panel_report/), [workers/](workers/) and [runtime/](runtime/) | Library and target processing, manuscript panel models and archived workflow components |
| [tools/](tools/) | Optional manifest import, pilot runs and archive audits |

Run optional tools from the checkout root, for example `python -m tools.luminary_pbc --help`. In the v0.2.0 release these utilities are top-level scripts.

The analysis uses a pinned revision of the upstream code and reference assets. [Upstream ABSOLVE-seq repository](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4) · [Upstream license](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/LICENSE.md)
