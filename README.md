<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/absolve-seq-agent-dark.png" />
    <img src="docs/assets/absolve-seq-agent.png" width="480" alt="ABSOLVE-seq Agent" />
  </picture>
</p>

<h1 align="center">From ABSOLVE-seq FASTQs to an indel report</h1>

<p align="center">
  Analyze the original 79-target reporter panel by talking to Codex or Claude Code.<br />
  Run on a workstation or Slurm. Open the report and follow each result back to its counts.
</p>

<p align="center">
  <a href="#start-with-one-message">💬 <strong>Start in chat</strong></a> &nbsp;·&nbsp;
  <a href="docs/d3734-t1-quick-start.md">🧬 Public dataset tutorial</a> &nbsp;·&nbsp;
  <a href="#example-conversations">📖 Example conversations</a> &nbsp;·&nbsp;
  <a href="#figure-3b-and-reproducibility">📊 Figure 3b &amp; Excel</a>
</p>

## Start with one message

**New to the agent?** Open an empty project folder in a signed-in Codex or Claude Code session on a machine with Git and internet access. Paste this message to set up the project **and run the public example**:

> Set up https://github.com/JasonLinjc/ABSOLVE-seq-agent in this empty folder. Follow the public D3734/T1 quick start, prepare the fast-mode environment and references, and download the paired NoEP and HiFi FASTQs. Run all 79 targets with the defaults, verify the results against the public benchmark, and give me the report and result table.

This example uses four compressed FASTQs (1.58 GB) from [SRA study SRP722226](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP722226). A completed run produces **316 target/strategy comparisons**. Open `results/report.md` for read retention and a four-strategy summary, or `results/paired_target_rates.tsv` for the counts behind any target.

**Already have original-study FASTQs?** Once the project is set up, give the agent the paths visible on the machine that will run the analysis:

> My FASTQs are in `/path/to/fastqs` and my sample sheet is at `/path/to/sample-sheet`. Find matched HiFi and NoEP samples. If more than one donor/replicate pair is available, ask me which to use. Run the selected pair with the study defaults and return the report, QC and target-level rates.

## Set up in Codex or Claude Code

Both skills ship with the repository. To set up manually, clone it and open the checkout in your client.

```bash
git clone https://github.com/JasonLinjc/ABSOLVE-seq-agent.git
cd ABSOLVE-seq-agent
```

| Client | Open the project | Try the analysis skill |
| --- | --- | --- |
| Codex | Open the cloned folder in Codex, or select Codex in the ChatGPT desktop app. For the CLI, run `codex` in the checkout and sign in with ChatGPT. | Type `$absolve-seq` in the CLI, or ask for the skill by name. |
| Claude Code | Run `claude` in the cloned folder and sign in. | Enter `/absolve-seq` followed by your request. |

The clients use the same [analysis instructions](.agents/skills/absolve-seq/SKILL.md) and [Figure 3b instructions](.agents/skills/absolve-seq-figures/SKILL.md). For client installation, see the [Codex quick start](https://developers.openai.com/codex/quickstart/) or [Claude Code quick start](https://code.claude.com/docs/en/quickstart).

<details>
<summary><strong>🔧 Check setup before downloading reads</strong></summary>

Ask: “List the pairs in `examples/manifest.example.tsv` and stop there.” You should see D3734/T1 with `NoEP_pool_D3734_T1` and `HiFi_pool_D3734_T1`. The example manifest contains placeholder paths.

```bash
ABSOLVE_MANIFEST=examples/manifest.example.tsv python paired.py --list
```

After installing the analysis environment, check that CRISPResso reports 2.3.4 and that the model and plotting dependencies load:

```bash
CRISPResso --version
Rscript -e 'stopifnot(requireNamespace("brglm2", quietly=TRUE))'
python -c 'import matplotlib, numpy, openpyxl'
```

For a real run, follow the [reference, manifest and FASTQ setup](docs/d3734-t1-quick-start.md). For the pinned Linux environment, use the [Docker or Apptainer instructions](docs/container.md).

</details>

## How the agent works

Tell the agent which samples you want to compare and where your FASTQs and sample sheet are, or ask it to use the public D3734/T1 example. It checks the read pairs and references, matches samples by donor and replicate, and asks you to choose if the comparison is ambiguous. You can select how the Maxi plasmid pool defines a qualifying barcode (pBC). The default follows the manuscript.

```mermaid
flowchart TD
    reads["FASTQs and sample sheet"] --> align["Prepare reads and call indels with CRISPResso2"]
    align --> counts["Count reads and pBC families under four strategies"]
    counts --> pair["Compare HiFi with its matched NoEP control"]
    counts --> panel["Fit WT and HiFi effects across donors"]
    frozen["Frozen library counts"] --> panel
    pair --> report["Verify the tables and deliver the report"]
    panel --> report
```

The agent runs locally or on Slurm, waits for the analysis to finish, and checks that every expected comparison appears once. For a recognized benchmark with matching inputs and settings, it also verifies the result tables against their recorded SHA-256 checksums.

| Your comparison | What to open |
| --- | --- |
| 🧬 **One donor and replicate** | Open `results/report.md` for QC and a summary. In `paired_target_rates.tsv`, each of the **316 rows** gives the HiFi and NoEP counts, both indel rates, and their difference in percentage points for one target and strategy. |
| 📋 **Full manuscript panel** | Supply all 26 libraries and submitted Supplementary Tables 2–5, or start from frozen `library_counts.tsv` and the submitted references. Open `report.md` in the results directory. Its `panel_comparison.tsv` covers **632 comparisons** across 79 targets, four strategies, and two editors, with adjusted q-values and agreement with the submitted calls. |

📊 **To reproduce Figure 3b**, ask for the submitted estimates or your panel estimates and provide Supplementary Table 2 for power. The agent returns the PDF and PNG with `Figure3b_source_data.xlsx`, which contains the 158 plotted target–editor rows. [Figure 3b inputs and outputs](docs/manuscript-figures.md)

## Example conversations

### From paired FASTQs to a report

The first prompt above runs the whole example. To guide it step by step, send these messages in the same chat:

1. **Identify inputs.** “Show me which mates belong to NoEP and HiFi for public D3734/T1. Verify their archive checksums and build my sample manifest.”
2. **Run the pair.** “Use four local workers, fast read processing and the manuscript Maxi pBC rule. Save all four strategies under `runs/D3734-T1` and wait for completion.”
3. **Check a result.** “Compare raw and prepared read pairs in the two arms. Does `1617_OT_0040_ALT` still have a positive difference after Dedud+Dedup? Show the counts.”
4. **Open the report.** “Give me the report and target table. Flag undefined rates and tell me whether the D3734/T1 benchmark passed.”

**What a verified answer looks like.** For that target in the public pair, NoEP has **0/1,619** edited pBC families and HiFi has **9/675**. The rates are 0.00% and 1.33%, giving **+1.33 percentage points** for HiFi minus NoEP. You can check both counts in `paired_target_rates.tsv`.

### From the full cohort to a manuscript report

> I have all 26 libraries, their manifest and submitted Supplementary Tables 2–5. Check the historical software versions, run the manuscript workflow and summarize agreement for WT and HiFi by strategy. Draw Figure 3b from the fitted Dedud estimates, use Supplementary Table 2 for power, and link the report, PDF and Excel source data.

If you already have frozen counts, ask the agent to start with a [count refit](docs/historical-maxi.md#refit-existing-counts) and identify that starting point in the report.

## Quick start

Prefer the terminal? Follow the [D3734/T1 walkthrough](docs/d3734-t1-quick-start.md) to install the environment, fetch the pinned upstream references, verify both mates and create the manifest. The public runs are NoEP `SRR39887646` and HiFi `SRR39887695`. The same donor and replicate also have a [WT library](https://www.ebi.ac.uk/ena/browser/view/SRR39887693) for the [full-panel comparison](#manuscript-comparison).

With inputs ready, run the public pair on a workstation:

```bash
python paired.py --donor D3734 --replicate T1 --threads 4 --fresh --output runs/D3734-T1
python verify_results.py pair-fast runs/D3734-T1/results
```

Local execution defaults to one worker. Use `--executor local` to select local execution on a host with Slurm tools installed. On a cluster, submit `sbatch paired.sbatch --donor D3734 --replicate T1` from the repository root in the analysis environment.

## Paired results

Each row gives the NoEP and HiFi numerators, denominators, rates (%) and difference (percentage points). A zero denominator leaves the rate blank. The fast workflow also writes `qc.json` and `strategy_summary.tsv` alongside the report.

| Strategy | What counts toward the rate |
| --- | --- |
| Unfiltered | All aligned reads |
| Dedup | One representative allele per sample, target and pBC family |
| Dedud | Reads with a qualifying Maxi pBC that pass cross-target filtering |
| Dedud+Dedup | One representative allele per retained family after filtering |

The default `--maxi-method manuscript` uses the published plasmid-barcode (pBC) annotation. The historical method rebuilds membership from six original Maxi FASTQs. [Historical Maxi reference and validation](docs/historical-maxi.md)

<details>
<summary><strong>Compare pBC definitions or change the target filter</strong></summary>

Select `merged` for QC on pooled Maxi reads, or `any-pool`, `two-pools` or `all-pools` to require support from one, two or all three pools. Set `ABSOLVE_MAXI_MERGED` to the merged reference or `ABSOLVE_MAXI_POOLS` to the parent of `LV_Maxi_1`, `LV_Maxi_2` and `LV_Maxi_3`. Each reference needs its `summary.json` and checksum-matched `plasmid_barcode_category.tsv.gz`.

The manuscript cross-target filter requires edited alleles to match the assigned target uniquely by minimum edit distance. Use `--target-filter pbc-only` for earlier barcode-only fast counts. To use the rebuilt reference, select `--maxi-method historical --maxi-reference DIR`. Run each comparison in a separate output directory.

</details>

## Figure 3b and reproducibility

Ask for **Figure 3b from the submitted workbook** to redraw the submitted estimates, or **Figure 3b from my panel refit** to use the agent's estimates. Both use candidate-specific power from Supplementary Table 2. The Excel workbook and plotted TSV contain one Dedud row per target and editor. [Figure commands and verification](docs/manuscript-figures.md)

The separate benchmark profiles cover a D3734/T1 fast pair, a manuscript-read pair with historical Maxi membership, and a frozen-count refit. Run the profile that matches your inputs and settings:

```bash
python verify_results.py pair-fast runs/paired-fast/results
python verify_results.py pair-historical runs/paired-historical/results
python verify_results.py panel-refit runs/panel-refit
```

## Manuscript comparison

The full design covers 26 libraries and 79 targets. D3734/T3 has WT and NoEP without a matched HiFi library. Manuscript read processing uses Trimmomatic 0.36 and paired-end CRISPResso2 2.2.13. [Historical runtime and input setup](docs/manuscript-reads.md)

```bash
python panel_report/references.py /path/to/submitted/workbooks data/submitted_model_reference --full-model
python paired.py --manuscript-panel --read-mode manuscript --output runs/manuscript-replay
```

The manuscript panel fits a Firth binomial regression with treatment and donor, then averages the three donor-specific logits equally for each treatment. A primary call requires Δ > 0.1 percentage points and BH q < 0.05, adjusted across 79 targets per editor and strategy. The fast panel weights donors by observed pBC-family counts.

| Validation starting point | Agreement with submitted primary calls |
| --- | --- |
| Historical-trimming FASTQ replay | 628/632 (four WT differences). [Count-level audit](docs/manuscript-discordance.md) |
| Frozen encounter-order counts, equal-donor refit | 632/632. [Refit validation and numerical differences](docs/historical-maxi.md#luminary-acceptance) |

## Local inputs and caches

Keep FASTQs, sample manifests, submitted workbooks and run outputs on storage accessible to the workflow. `ABSOLVE_DATA` selects the data directory, `ABSOLVE_MANIFEST` selects a manifest and `ABSOLVE_CACHE` points to compatible cached alignments. `--fresh` forces a paired run to realign. Use the [manifest template](examples/manifest.example.tsv) for laboratory files, or build it with `python resolve_inputs.py INVENTORY.tsv METADATA.xlsx data/sample_manifest.tsv`.

The analysis uses a pinned revision of the [upstream ABSOLVE-seq repository](https://github.com/pellinlab/ABSOLVE-seq/tree/271355dfe40cf46b952e82b88dd4521f23f9feb4). Its code and reference assets are covered by the [upstream license](https://github.com/pellinlab/ABSOLVE-seq/blob/271355dfe40cf46b952e82b88dd4521f23f9feb4/LICENSE.md).
