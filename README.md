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

**Try the public example.** Open an empty folder in Codex or Claude Code and paste this message.

> Set up https://github.com/JasonLinjc/ABSOLVE-seq-agent and run the public D3734/T1 example from the raw sequencing reads. How much editing do we see with HiFi compared with its control at each target? Show me the report, figure and Excel data.

The agent downloads four FASTQ files (1.58 GB) and analyzes all 79 targets using both ends of each read pair and the published plasmid barcode annotation. Before handing you the report, it checks the editing rates against the verified example and checks the power estimates for detecting small editing differences. Start with `results/report.md` for the figure and links to the tables. [Public dataset](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP722226) · [Verified example results](docs/d3734-t1-quick-start.md#verified-on-luminary) · [Docker acceptance](docs/container.md#d3734t1-docker-acceptance)

You’ll need a signed-in client, Git, internet access and either Conda/Mamba or a supported container runtime. [Download the Docker or Apptainer image](https://github.com/JasonLinjc/ABSOLVE-seq-agent/releases/tag/v0.2.0) · [Container setup](docs/container.md#download-a-prebuilt-image)

**Using a cluster from your laptop?** Add the cluster’s SSH name and an empty shared folder to your message so the agent knows where to run the analysis. [Cluster setup and transfer](docs/container.md#stage-public-inputs-from-a-networked-machine)

**Already have original-study FASTQs?** Once the project is set up, give the agent the paths visible on the machine that will run the analysis:

> My sequencing files are in `/path/to/fastqs` and my sample sheet is at `/path/to/sample-sheet`. Find a HiFi sample and its NoEP control from the same donor and replicate. If there are several pairs, let me choose one. Then analyze it and show me the results.

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

In the paired analysis environment, check that CRISPResso reports 2.2.13 and that the model and plotting dependencies load.

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

The default workflow trims both mates with Trimmomatic 0.36 and assigns each read pair to a target by tBC. **CRISPResso2 2.2.13** uses FLASH to merge the mates, aligns the merged reads to each target’s reference amplicon and identifies insertions and deletions in the cut-site window. The agent then groups aligned reads by pBC and calculates indel rates under the four counting strategies.

<p align="center">
  <img src="docs/assets/agent-architecture.svg" width="660" alt="A chat request selects the analysis or figure skill in Codex or Claude Code. The skills use the CLI directly. The published pBC annotation is the default. The original study’s six Maxi FASTQs can rebuild the annotation under the published or historical rule for reuse across runs. The local or Slurm workflow assigns reads by tBC, uses CRISPResso2 2.2.13 to merge and align the mates and identify cut-site indels, groups them by pBC, estimates NoEP and HiFi indel rates, calculates Δindel rate (HiFi minus NoEP) and simulates detection power. It exports indel rates and paired differences in paired_target_rates.tsv and detection power in power_estimates.tsv, alongside a report and a figure with Excel source data." />
</p>

**What to prepare.** Provide four FASTQ files (R1 and R2 for both NoEP and HiFi) from the same donor and replicate, plus a sample sheet or manifest that maps each file to its sample. The agent uses the published pBC annotation by default. To rebuild it from reads, use the original study’s six Maxi FASTQs (R1 and R2 for each of three pools) and choose the published or historical rule. Later runs can reuse the saved annotation directory. To try the agent with public data, ask for the [D3734/T1 example](docs/d3734-t1-quick-start.md).

**Where to find results.** Open `results/report.md` in your run directory. In the same folder, `paired_target_rates.tsv` gives indel rates and Δindel rates (HiFi minus NoEP), and `power_estimates.tsv` gives detection power by target and strategy. The PDF, PNG and Excel source data are in `results/figures/`.

## Example conversations

### From paired FASTQs to a report

The [one-message prompt](#start-with-one-message) runs the public example end to end. Here is the same analysis step by step, with results from the [verified D3734/T1 run](docs/d3734-t1-quick-start.md#verified-on-luminary).

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Can you get the public D3734/T1 example? Which sample is the control?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>NoEP is the control (<code>SRR39887646</code>). Its matched HiFi sample is <code>SRR39887695</code>. I checked all four sequencing files against the archive and recorded them in the sample list.</td></tr>
<tr><td><strong>👤 User</strong><br>Analyze all 79 targets from those reads and show me how the two samples compare.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>I used both ends of each read pair and the published pBC annotation to compare HiFi with its matched NoEP control at all 79 targets. The report shows 316 comparisons, one per target and counting strategy.</td></tr>
<tr><td><strong>👤 User</strong><br>At <code>1617_OT_0040_ALT</code>, do we still see more editing with HiFi if we count each qualifying plasmid barcode just once?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Yes. With that counting method (Dedud+Dedup), HiFi had 11 edited barcode families out of 691 (1.59%). NoEP had 0 out of 1,694 (0.00%). The difference was <strong>+1.59 percentage points</strong>.</td></tr>
<tr><td><strong>👤 User</strong><br>Did my results match the public example? Where can I see the figure and the numbers behind it?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>Yes. All 316 editing comparisons matched the verified example. I also checked the power estimates and the data used for the figure. Open <code>runs/D3734-T1/results/report.md</code> for the figure and links to the editing rates, power estimates and Excel file.</td></tr>
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

<details>
<summary><strong>Reproduce Supplementary Tables 2–5</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>I have the original 26 libraries and the submitted workbooks. Can you reproduce Supplementary Tables 2–5 and show me any differences?</td></tr>
<tr><td><strong>🤖 Agent</strong><br>I’ll use the historical processing recorded for WT, NoEP and HiFi. I’ll check the historical pBC annotation, compare all 632 target, strategy and editor combinations with the workbooks, and report differences in rates, q-values and positive calls.</td></tr>
</tbody>
</table>

</details>

<a id="validate-pbc-annotation-reproducibility"></a>
<details>
<summary><strong>Validate pBC annotation reproducibility</strong></summary>

<table>
<tbody>
<tr><td><strong>👤 User</strong><br>Validate pBC annotation reproducibility.</td></tr>
<tr><td><strong>🤖 Agent</strong><br>I’ll compare the published and rebuilt annotations against the saved alignments for all 26 libraries in the Luminary archive. I’ll return the comparison summary, field-level differences and SHA-256 manifest.</td></tr>
</tbody>
</table>

You can make this request in your preferred language in Codex or Claude Code. [One-time setup and output guide](docs/pbc-validation.md#run-from-codex-or-claude-code)

</details>

## Quick start

Prefer the terminal? Follow the [D3734/T1 walkthrough](docs/d3734-t1-quick-start.md) to install the environment, fetch the pinned upstream references, verify both mates and create the manifest. The public runs are NoEP `SRR39887646` and HiFi `SRR39887695`. The same donor and replicate also have a [WT library](https://www.ebi.ac.uk/ena/browser/view/SRR39887693) for the [full-panel comparison](#manuscript-comparison).

With inputs ready, run the public pair on a workstation:

```bash
python paired.py --executor local --donor D3734 --replicate T1 --threads 4 --output runs/D3734-T1
python verify_results.py pair-manuscript runs/D3734-T1/results
python verify_results.py power runs/D3734-T1/results
```

The first command starts from the raw reads and writes a new run directory. The next two commands check the public example against its 316 accepted results and verify that the power estimates, plot data and Excel values agree.

On Slurm, run `python paired.py --donor D3734 --replicate T1 --output runs/D3734-T1` on the submission host with the analysis environment active. It prints the preparation, alignment and finish job IDs. Wait for the finish job before opening the report. [Cluster setup and verification](docs/d3734-t1-quick-start.md#slurm-submission)

## Paired results

In `results/paired_target_rates.tsv`, each row gives the NoEP and HiFi numerators, denominators, rates (%) and difference (percentage points). A zero denominator leaves the rate blank. The default workflow writes `library_counts.tsv`, `sample_inventory.tsv` and `provenance.json` alongside the report.

| Strategy | What counts toward the rate |
| --- | --- |
| Unfiltered | All aligned reads |
| Dedup | One most frequent allele per sample, target and pBC |
| Dedud | Reads with a qualifying plasmid barcode that pass cross-target filtering |
| Dedud+Dedup | One representative allele per retained family after filtering |

When alleles tie for the highest read count, the paired-end workflow breaks the tie by aligned read sequence and then aligned reference sequence.

Plasmid pool sequencing identifies which plasmid barcodes (pBCs) qualify for each target. The default `--maxi-method manuscript` uses the published annotation derived from the study’s Maxi plasmid pool. To rebuild from the six original Maxi FASTQs, use `setup_references.py --published-maxi` for the published entropy threshold or `--historical-maxi` for the historical two-tBC exception. On later runs, pass the saved directory with `--maxi-reference DIR` and its matching `--maxi-method manuscript` or `--maxi-method historical`. [Build commands and rule comparison](docs/historical-maxi.md#build-the-reference-once)

[Full-panel pBC annotation validation using archived alignments](docs/pbc-validation.md)

## Figure 3b and reproducibility

The panel workflow writes `figures/Figure3b.pdf` and `figures/Figure3b_source_data.xlsx` using its fitted effects and freshly simulated power. To reproduce the published graphic, ask for **Figure 3b with submitted power**. The redraw can use submitted or refitted effects while retaining the original Supplementary Table 2 power values. Both versions export 158 plotted target/editor rows and record the power source. [Figure commands and verification](docs/manuscript-figures.md)

## Manuscript comparison

To reproduce the paper, ask “Reproduce Supplementary Tables 2–5 from the original 26 libraries.” The full-panel command `python paired.py --manuscript-panel` selects `historical` processing for each study arm. Supply the historical pBC annotation and WT runtime as described in the guide below.

[Raw FASTQ reproduction of Supplementary Tables 2–5: acceptance results and run instructions](docs/historical-panel-acceptance.md)

## Local inputs and caches

Keep inputs on storage visible to the workers. `ABSOLVE_DATA` selects the data directory (default `data/` in the checkout). `ABSOLVE_MANIFEST` defaults to `$ABSOLVE_DATA/sample_manifest.tsv`. Use a new output directory for each paired-end run.

For laboratory files, follow the [manifest template](examples/manifest.example.tsv). Supply the original study sample names, absolute R1 and R2 paths and their SHA-256 digests. The importer `resolve_inputs.py INVENTORY.tsv METADATA.xlsx OUTPUT.tsv` joins a checksum-verified inventory to the submission workbook’s `SRA_data` sheet by both original filenames. The public quick start builds a manifest for the D3734/T1 NoEP and HiFi pair from the four ENA downloads.

The analysis uses a pinned revision of the [upstream ABSOLVE-seq repository](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4). Its code and reference assets are covered by the [upstream license](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/LICENSE.md).
