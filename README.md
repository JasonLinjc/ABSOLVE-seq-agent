<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/absolve-seq-agent-dark.png" />
    <img src="docs/assets/absolve-seq-agent.png" width="480" alt="ABSOLVE-seq Agent" />
  </picture>
</p>

<h1 align="center">Quantify off-target gene editing with ABSOLVE-seq</h1>

<p align="center">
  Ask Codex or Claude Code to analyze your FASTQs from the original 79-target reporter panel.<br />
  Get indel rates, paired editing effects and detection power in a report with a figure and Excel source data.
</p>

<p align="center">
  <a href="#start-with-one-message">💬 <strong>Start in chat</strong></a> &nbsp;·&nbsp;
  <a href="docs/d3734-t1-quick-start.md">🧬 Public dataset tutorial</a> &nbsp;·&nbsp;
  <a href="#example-conversations">📖 Example conversations</a> &nbsp;·&nbsp;
  <a href="#figure-3b-and-reproducibility">📊 Figure 3b &amp; Excel</a>
</p>

**What is ABSOLVE-seq?** ABSOLVE-seq (Assessment By Stand-in Off-target LentiViral Ensemble with sequencing) measures how often CRISPR causes insertions or deletions (indels) in reporters carrying reference or variant sequences from candidate off-target sites. From FASTQs, the agent reports indel rates for HiFi Cas9 and its matched control, their difference at each target, and estimated detection power, with a figure and Excel table. [Original ABSOLVE-seq study](https://doi.org/10.1101/2024.07.24.605019)

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

<a id="fastq-to-report-workflow"></a>

### Agent architecture

After target assignment by tBC, **CRISPResso2** aligns reads to each target’s reference amplicon and identifies insertions and deletions in the cut-site window. The agent then groups aligned reads by pBC and calculates indel rates under the four counting strategies.

<p align="center">
  <img src="docs/assets/agent-architecture.svg" width="660" alt="A chat request selects the analysis or figure skill in Codex or Claude Code. The skills use the CLI directly. The published pBC annotation is the default. The original study’s six Maxi FASTQs can rebuild the annotation under the published or historical rule for reuse across runs. The local or Slurm workflow assigns reads by tBC, uses CRISPResso2 to align them and identify cut-site indels, groups them by pBC, estimates NoEP and HiFi indel rates, calculates Δindel rate (HiFi minus NoEP) and simulates detection power. It exports indel rates and paired differences in paired_target_rates.tsv and detection power in power_estimates.tsv, alongside a report and a figure with Excel source data." />
</p>

**What to prepare.** Provide four FASTQ files (R1 and R2 for both NoEP and HiFi) from the same donor and replicate, plus a sample sheet or manifest that maps each file to its sample. The agent uses the published pBC annotation by default. To rebuild it from reads, use the original study’s six Maxi FASTQs (R1 and R2 for each of three pools) and choose the published or historical rule. Later runs can reuse the saved annotation directory. To try the agent with public data, ask for the [D3734/T1 example](docs/d3734-t1-quick-start.md).

**Where to find results.** Open `results/report.md` in your run directory. In the same folder, `paired_target_rates.tsv` gives indel rates and Δindel rates (HiFi minus NoEP), and `power_estimates.tsv` gives detection power by target and strategy. The PDF, PNG and Excel source data are in `results/figures/`.

## Example conversations

### From paired FASTQs to a report

The [one-message prompt](#start-with-one-message) runs the public example end to end. Here is the same analysis step by step, with results from the [verified D3734/T1 run](docs/d3734-t1-quick-start.md#verified-on-luminary).

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Find the paired FASTQs for public D3734/T1. Verify their archive checksums and build my sample manifest.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>NoEP is <code>SRR39887646</code> and HiFi is <code>SRR39887695</code>. Both mates of each run passed the archive checksum check. The manifest records all four FASTQs and their SHA-256 digests.</td></tr>
<tr><td><strong>👤 User</strong><br>Run all 79 targets from scratch with four local workers, fast read processing and the published pBC annotation. Save the results under <code>runs/D3734-T1</code> and wait for completion.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>The run produced 316 comparisons across four strategies. Read preparation retained 5,971,199 of 10,084,551 NoEP pairs and 6,979,407 of 9,868,448 HiFi pairs.</td></tr>
<tr><td><strong>👤 User</strong><br>Is the HiFi minus NoEP difference positive at <code>1617_OT_0040_ALT</code> after Dedud+Dedup? Show me the counts.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Yes. HiFi had 9 edited pBC families out of 675 (1.33%). NoEP had 0 out of 1,619 (0.00%). The difference was <strong>+1.33 percentage points</strong>.</td></tr>
<tr><td><strong>👤 User</strong><br>Check the pair benchmark and power estimates. Where are my report, tables, figure and Excel data?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Both checks passed. Open <code>runs/D3734-T1/results/report.md</code>. In the same folder, <code>paired_target_rates.tsv</code> contains the target-level rates and differences, and <code>power_estimates.tsv</code> contains the power estimates. The figure and its Excel source data are in <code>figures/</code>.</td></tr>
</tbody>
</table>

<a id="choose-settings-in-chat"></a>
<details>
<summary><strong>Choose a pBC annotation</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Use the public D3734/T1 pair on my workstation. Let me choose the pBC annotation before you start.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Which pBC annotation should I use for Dedud and Dedud+Dedup?
<ol>
<li>Use the published annotation (default).</li>
<li>Rebuild from the original study’s six Maxi FASTQs. Use the published rule, or select the historical two-tBC exception.</li>
<li>Reuse the annotation from a previous Maxi rebuild with its saved rule.</li>
</ol>
</td></tr>
<tr><td><strong>👤 User</strong><br>Option 3. I rebuilt it on my last run.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Where did you save the annotation directory?</td></tr>
<tr><td><strong>👤 User</strong><br>Use the annotation in <code>/data/annotations/maxi</code> and run the pair.</td></tr>
</tbody>
</table>

You can also ask to choose among donor/replicate pairs in your manifest or to run on Slurm. For option 2, the agent checks the six compressed FASTQs against the [study file manifest](runtime/maxi/inputs.sha256).

</details>

<a id="validate-pbc-annotation-reproducibility"></a>
<details>
<summary><strong>Validate pBC annotation reproducibility</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Validate pBC annotation reproducibility.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>I’ll use the configured Luminary archive to compare the published and rebuilt annotations across all 26 libraries in fast and manuscript reproduction modes. I’ll return the comparison summary, field-level differences and SHA-256 manifest.</td></tr>
</tbody>
</table>

You can make this request in your preferred language in Codex or Claude Code. [One-time setup and output guide](docs/pbc-validation.md#run-from-codex-or-claude-code)

</details>

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

Plasmid pool sequencing identifies which plasmid barcodes (pBCs) qualify for each target. The default `--maxi-method manuscript` uses the published annotation derived from the study’s Maxi plasmid pool. To rebuild from the six original Maxi FASTQs, use `setup_references.py --published-maxi` for the published entropy threshold or `--historical-maxi` for the historical two-tBC exception. On later runs, pass the saved directory with `--maxi-reference DIR` and its matching `--maxi-method manuscript` or `--maxi-method historical`. [Build commands and rule comparison](docs/historical-maxi.md#build-the-reference-once)

[Full-panel pBC validation in fast and manuscript reproduction modes](docs/pbc-validation.md)

<details>
<summary><strong>Compare pBC definitions or change the target filter</strong></summary>

For a **fast matched pair**, select `merged` for QC on pooled Maxi reads, or `any-pool`, `two-pools` or `all-pools` to require support from one, two or all three pools. Set `ABSOLVE_MAXI_MERGED` to the merged reference or `ABSOLVE_MAXI_POOLS` to the parent of `LV_Maxi_1`, `LV_Maxi_2` and `LV_Maxi_3`. Each reference needs its `summary.json` and checksum-matched `plasmid_barcode_category.tsv.gz`.

The manuscript cross-target filter requires edited alleles to match the assigned target uniquely by minimum edit distance. A fast pair can use `--target-filter pbc-only` to skip this check. Manuscript-read pairs and both full-panel modes always apply the cross-target check and accept the published annotation or a Maxi FASTQ rebuild made with either rule.

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

[Raw FASTQ reproduction of Supplementary Tables 2–5: acceptance results and run instructions](docs/historical-panel-acceptance.md)

## Local inputs and caches

Keep inputs on storage visible to the workers. `ABSOLVE_DATA` selects the data directory (default `data/` in the checkout). `ABSOLVE_MANIFEST` defaults to `$ABSOLVE_DATA/sample_manifest.tsv`. A fast pair can reuse a compatible alignment registered in the study cache under `ABSOLVE_CACHE`. Use `--fresh` to prepare and align the FASTQs again.

For laboratory files, follow the [manifest template](examples/manifest.example.tsv). Supply the original study sample names, absolute R1 and R2 paths and their SHA-256 digests. The importer `resolve_inputs.py INVENTORY.tsv METADATA.xlsx OUTPUT.tsv` joins a checksum-verified inventory to the submission workbook’s `SRA_data` sheet by both original filenames. The public quick start builds a manifest for the D3734/T1 NoEP and HiFi pair from the four ENA downloads.

The analysis uses a pinned revision of the [upstream ABSOLVE-seq repository](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4). Its code and reference assets are covered by the [upstream license](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/LICENSE.md).
