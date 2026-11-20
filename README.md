<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/absolve-seq-agent-dark.png" />
    <img src="docs/assets/absolve-seq-agent.png" width="480" alt="ABSOLVE-seq Agent" />
  </picture>
</p>

<h1 align="center">From ABSOLVE-seq FASTQs to an indel report</h1>

<p align="center">
  Analyze the original 79-target reporter panel by talking to Codex or Claude Code.<br />
  Run on a workstation or Slurm. Get editing effects, simulated detection power, a report and a figure.
</p>

<p align="center">
  <a href="#start-with-one-message">💬 <strong>Start in chat</strong></a> &nbsp;·&nbsp;
  <a href="docs/d3734-t1-quick-start.md">🧬 Public dataset tutorial</a> &nbsp;·&nbsp;
  <a href="#example-conversations">📖 Example conversations</a> &nbsp;·&nbsp;
  <a href="#figure-3b-and-reproducibility">📊 Figure 3b &amp; Excel</a>
</p>

## Start with one message

**New to the agent?** Open an empty project folder in a signed-in Codex or Claude Code session on a machine with Git, internet access and either Conda/Mamba or a supported container runtime. Paste this message to set up the project **and run the public example**:

> Set up https://github.com/JasonLinjc/ABSOLVE-seq-agent in this empty folder. Follow the public D3734/T1 quick start, prepare the fast-mode environment and references, and download the paired NoEP and HiFi FASTQs. Run all 79 targets with `--fresh` and the default filters. Wait for completion, run the `pair-fast` and `power` checks in the analysis environment, and give me the report, effect and power tables, figure and Excel source data.

This example uses four compressed FASTQs (1.58 GB) from [SRA study SRP722226](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP722226). A completed run produces **316 target/strategy comparisons**. Open `results/report.md` for read retention and a four-strategy summary, or `results/paired_target_rates.tsv` for the counts behind any target. [Clean-checkout acceptance results](docs/d3734-t1-quick-start.md#verified-on-luminary)

**Using a cluster from your laptop?** Include its SSH host and an empty shared project directory in the same message. The agent can use a login shell to find Slurm and stage public downloads from a networked machine if cluster downloads are slow. It builds the reference index and verifies the four FASTQs in the new checkout before submitting the pair. [Cluster setup and transfer](docs/container.md#stage-public-inputs-from-a-networked-machine)

**Already have original-study FASTQs?** Once the project is set up, give the agent the paths visible on the machine that will run the analysis:

> My FASTQs are in `/path/to/fastqs` and my sample sheet is at `/path/to/sample-sheet`. Find matched HiFi and NoEP samples. If more than one donor/replicate pair is available, ask me which to use. Run the selected pair with the study defaults and return the report, target-level rates, estimated power and figure.

## Set up in Codex or Claude Code

Both skills ship with the repository. To set up manually, clone it and open the checkout in your client.

```bash
git clone https://github.com/JasonLinjc/ABSOLVE-seq-agent.git
cd ABSOLVE-seq-agent
```

| Client | Open the project | Try the analysis skill |
| --- | --- | --- |
| Codex | Open the cloned folder in the Codex app, or run `codex` in the checkout. Sign in with your ChatGPT account. | Type `$absolve-seq` in the CLI, or ask for the skill by name. |
| Claude Code | Run `claude` in the cloned folder and sign in. | Enter `/absolve-seq` followed by your request. |

The clients use the same [analysis instructions](.agents/skills/absolve-seq/SKILL.md) and [Figure 3b instructions](.agents/skills/absolve-seq-figures/SKILL.md). For client installation, see the [Codex quick start](https://developers.openai.com/codex/quickstart/) or [Claude Code quick start](https://code.claude.com/docs/en/quickstart).

<details>
<summary><strong>🔧 Check setup before downloading reads</strong></summary>

Ask: “List the pairs in `examples/manifest.example.tsv` and stop there.” You should see D3734/T1 with `NoEP_pool_D3734_T1` and `HiFi_pool_D3734_T1`. This command lists matched manifest entries and flags a missing arm. The example paths are placeholders. File and checksum checks happen during setup and processing.

```bash
ABSOLVE_MANIFEST=examples/manifest.example.tsv python paired.py --list
```

After installing the fast-mode environment, check that CRISPResso reports 2.3.4 and that the model and plotting dependencies load:

```bash
CRISPResso --version
Rscript -e 'stopifnot(requireNamespace("brglm2", quietly=TRUE))'
python -c 'import matplotlib, numpy, openpyxl'
```

For a real run, follow the [reference, manifest and FASTQ setup](docs/d3734-t1-quick-start.md). For the pinned Linux environment, use the [Docker or Apptainer instructions](docs/container.md).

</details>

## How the agent works

### Agent architecture

<p align="center">
  <img src="docs/assets/agent-architecture.svg" width="660" alt="A chat request selects the analysis or figure skill in Codex or Claude Code. The skills use the CLI directly. The software environment and sequencing inputs feed a local or Slurm workflow that returns a report, effect and power tables, and a figure with Excel source data." />
</p>

Give the agent the FASTQ folder and sample sheet for your NoEP and HiFi libraries, or ask for the public D3734/T1 example. It matches samples by donor and replicate and asks which pair to use if several match.

The analysis skill runs `paired.py` locally or on Slurm. The public setup indexes the target barcodes, reference amplicons, guides and plasmid barcode annotations used to process the reads. Codex and Claude Code call the CLI directly. [Analysis environment](docs/container.md)

**From reads to counts.** A matched pair defaults to fast read processing. The agent uses the **target barcode (tBC)** to assign reads to a target before alignment. Within each sample and target, reads with the same **plasmid barcode (pBC)** form one family. Unfiltered counts reads, Dedup counts one representative allele per family, and Dedud applies the selected plasmid barcode rule and, by default, a cross-target identity check. Dedud+Dedup combines filtering and family counting. The published plasmid barcode annotation from the Maxi pool is the default in both read modes. [Counting strategies](#paired-results) · [Read modes](docs/manuscript-reads.md)

**Editing effect.** For each target and strategy, the agent divides the indel count by the denominator in each sample, then subtracts the NoEP rate from the HiFi rate. Open `results/paired_target_rates.tsv` for the counts, rates and difference in percentage points. The table has **316 rows** (79 targets × 4 strategies).

**Detection power.** Using the observed depths and estimated NoEP background, the agent simulates 1,000 experiments at each of two assumed increases, **0.01 and 0.1 percentage points**. Power is the fraction of simulated one-sided tests with nominal *p* < 0.05. Open `results/power_estimates.tsv` for the 632 target/strategy/effect combinations, their status and Monte Carlo intervals. [Simulation method](docs/power.md)

**Report and figure.** Open `results/report.md` for the four-strategy summary and links to the Dedud figure (PDF and PNG) and its Excel source data. The analysis skill waits for completion and checks the outputs before returning them.

## Example conversations

### From paired FASTQs to a report

The first prompt above runs the whole example. To guide it step by step, send these messages in the same chat:

1. **Identify inputs.** “Show me which mates belong to NoEP and HiFi for public D3734/T1. Verify their archive checksums and build my sample manifest.”
2. **Run the pair.** “Use four local workers, fast read processing and the manuscript plasmid barcode rule. Save all four strategies under `runs/D3734-T1` and wait for completion.”
3. **Check a result.** “Compare raw and prepared read pairs in the two arms. Does `1617_OT_0040_ALT` still have a positive difference after Dedud+Dedup? Show the counts.”
4. **Open the report.** “Give me the report, figure, target table and power estimates. Flag undefined rates and tell me whether the D3734/T1 benchmark passed.”

**What a verified answer looks like.** For that target in the public pair, NoEP has **0/1,619** edited pBC families and HiFi has **9/675**. The rates are 0.00% and 1.33%, giving **+1.33 percentage points** for HiFi minus NoEP. You can check both counts in `paired_target_rates.tsv`.

## Quick start

Prefer the terminal? Follow the [D3734/T1 walkthrough](docs/d3734-t1-quick-start.md) to install the environment, fetch the pinned upstream references, verify both mates and create the manifest. The public runs are NoEP `SRR39887646` and HiFi `SRR39887695`. The same donor and replicate also have a [WT library](https://www.ebi.ac.uk/ena/browser/view/SRR39887693) for the [full-panel comparison](#manuscript-comparison).

With inputs ready, run the public pair on a workstation:

```bash
python paired.py --executor local --donor D3734 --replicate T1 --threads 4 --fresh --output runs/D3734-T1
python verify_results.py pair-fast runs/D3734-T1/results
python verify_results.py power runs/D3734-T1/results
```

Use a new output directory for each run. Local execution defaults to one worker. On Slurm, submit `sbatch paired.sbatch --donor D3734 --replicate T1 --fresh --output runs/D3734-T1` from the repository root. Run the two verification commands after the job completes.

## Paired results

In `results/paired_target_rates.tsv`, each row gives the NoEP and HiFi numerators, denominators, rates (%) and difference (percentage points). A zero denominator leaves the rate blank. Fast pairs also write `qc.json` and `strategy_summary.tsv`. Manuscript-read pairs write `library_counts.tsv`, `sample_inventory.tsv` and `provenance.json` alongside the report.

| Strategy | What counts toward the rate |
| --- | --- |
| Unfiltered | All aligned reads |
| Dedup | One most frequent allele per sample, target and pBC |
| Dedud | Reads with a qualifying plasmid barcode that pass cross-target filtering |
| Dedud+Dedup | One representative allele per retained family after filtering |

When alleles tie for the highest read count, fast reads use the first allele encountered. Manuscript reads break the tie by aligned read sequence and then aligned reference sequence.

The default `--maxi-method manuscript` uses the published plasmid barcode (pBC) annotation. To use historical membership, first build the reference from the six original Maxi FASTQs with `python setup_references.py --historical-maxi FASTQ_DIR --trimmomatic-jar JAR --output DIR`, then pass `--maxi-method historical --maxi-reference DIR` to the analysis. [Historical Maxi reference and validation](docs/historical-maxi.md)

<details>
<summary><strong>Compare pBC definitions or change the target filter</strong></summary>

For a **fast matched pair**, select `merged` for QC on pooled Maxi reads, or `any-pool`, `two-pools` or `all-pools` to require support from one, two or all three pools. Set `ABSOLVE_MAXI_MERGED` to the merged reference or `ABSOLVE_MAXI_POOLS` to the parent of `LV_Maxi_1`, `LV_Maxi_2` and `LV_Maxi_3`. Each reference needs its `summary.json` and checksum-matched `plasmid_barcode_category.tsv.gz`.

The manuscript cross-target filter requires edited alleles to match the assigned target uniquely by minimum edit distance. A fast pair can use `--target-filter pbc-only` to skip this check. Manuscript-read pairs and both full-panel modes always apply the cross-target check and accept either published (`manuscript`) or rebuilt (`historical`) pBC membership.

</details>

## Figure 3b and reproducibility

The panel workflow writes `figures/Figure3b.pdf` and `figures/Figure3b_source_data.xlsx` using its fitted effects and freshly simulated power. To reproduce the published graphic, ask for **Figure 3b with submitted power**. The redraw can use submitted or refitted effects while retaining the original Supplementary Table 2 power values. Both versions export 158 plotted target/editor rows and record the power source. [Figure commands and verification](docs/manuscript-figures.md)

The separate benchmark profiles cover a D3734/T1 fast pair, a manuscript-read pair with historical Maxi membership, and a panel refit from saved counts. Run the profile that matches your inputs and settings:

```bash
python verify_results.py pair-fast runs/paired-fast/results
python verify_results.py pair-historical runs/paired-historical/results
python verify_results.py panel-refit runs/panel-refit
```

## Manuscript comparison

The optional `--manuscript-panel` mode checks the agent against the submitted study results.

<details>
<summary><strong>Show reproduction inputs, commands and validation results</strong></summary>

> I have all 26 libraries, their manifest and submitted Supplementary Tables 2–5. Check the historical software versions, run the manuscript workflow and summarize agreement for WT and HiFi by strategy. Estimate power from my counts, draw Figure 3b from the fitted Dedud estimates and simulated power, and return the report, PDF and Excel source data.

If you already have a complete `library_counts.tsv` from a previous panel run, ask the agent to start with a [count refit](docs/historical-maxi.md#refit-existing-counts) and identify that starting point in the report.

The manuscript panel pools replicate counts within each donor and treatment, then fits a Firth logistic regression for each target and strategy. It averages the three donor-specific logits equally and converts the result to a rate for each treatment. A primary positive call requires an editor-minus-NoEP difference above **0.1 percentage points** and **q < 0.05**. The q-values come from one-sided Wald tests with Benjamini–Hochberg correction across 79 targets within each editor and strategy.

<p align="center">
  <img src="docs/assets/paired-workflow.svg" width="660" alt="The agent selects a matched pair or the full manuscript panel, checks inputs, assigns reads to targets by target barcode (tBC), aligns them with CRISPResso2, groups reads by plasmid barcode (pBC) and counts indels four ways. It estimates 316 paired differences or 632 panel contrasts, simulates power, exports the Dedud figure and Excel data, then checks the outputs and returns the report." />
</p>

The full design covers 26 libraries and 79 targets. D3734/T3 has WT and NoEP without a matched HiFi library. Manuscript read processing uses Trimmomatic 0.36 and paired-end CRISPResso2 2.2.13. [Historical runtime and input setup](docs/manuscript-reads.md)

```bash
python panel_report/references.py /path/to/submitted/workbooks "${ABSOLVE_DATA:-data}/submitted_model_reference" --full-model
python paired.py --manuscript-panel --read-mode manuscript --output runs/manuscript-replay
```

The first command reads `SupplementaryTable2.xlsx` through `SupplementaryTable5.xlsx`. Run it once in a fresh reference directory. On a workstation the second command runs to completion. On a Slurm submission host it records the **prepare**, **align** and **finish** job IDs in `runs/manuscript-replay/jobs.json`. Wait for **finish** to succeed before opening the report.

| Starting point | Effect estimate | Report path relative to `--output` |
| --- | --- | --- |
| Matched pair, fast or manuscript reads | Observed HiFi and NoEP rates for one donor and replicate | `results/report.md` |
| Full panel, manuscript reads (default) | Three donors weighted equally on the logit scale | `results/report.md` |
| Full panel, `--read-mode fast` | Donors weighted by observed pBC-family counts | `report.md` |
| Full panel, `--refit-counts FILE --output DIR` | Three donors weighted equally on the logit scale, using the filtering already present in the counts | `report.md` |

Open the tables linked from each report and the PDF, PNG and Excel files in the adjacent `figures/` directory. A count refit runs directly in the calling process. See the [refit command and required inputs](docs/historical-maxi.md#refit-existing-counts).

| Validation starting point | Agreement with submitted primary calls |
| --- | --- |
| Historical-trimming FASTQ replay | 628/632 (four WT differences). [Historical replay and four discordant calls](docs/manuscript-reads.md#luminary-validation) |
| Saved encounter-order counts, equal-donor refit | 632/632. [Refit validation and numerical differences](docs/historical-maxi.md#luminary-acceptance) |

</details>

## Local inputs and caches

Keep inputs on storage visible to the workers. `ABSOLVE_DATA` selects the data directory (default `data/` in the checkout). `ABSOLVE_MANIFEST` defaults to `$ABSOLVE_DATA/sample_manifest.tsv`. A fast pair can reuse a compatible alignment registered in the study cache under `ABSOLVE_CACHE`. Use `--fresh` to prepare and align the FASTQs again.

For laboratory files, follow the [manifest template](examples/manifest.example.tsv). Supply the original study sample names, absolute R1 and R2 paths and their SHA-256 digests. The importer `resolve_inputs.py INVENTORY.tsv METADATA.xlsx OUTPUT.tsv` joins a checksum-verified inventory to the submission workbook’s `SRA_data` sheet by both original filenames. The public quick start builds a manifest for the D3734/T1 NoEP and HiFi pair from the four ENA downloads.

The analysis uses a pinned revision of the [upstream ABSOLVE-seq repository](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4). Its code and reference assets are covered by the [upstream license](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/LICENSE.md).
