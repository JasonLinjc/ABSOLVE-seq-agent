# Acceptance of the historical panel from raw FASTQs

The Luminary run reproduced **all 632 classifications in Supplementary Tables 2–5 (63 positive and 569 negative)** from 52 checksum-verified raw FASTQs. It processed 334,970,765 read pairs from 26 libraries, including all nine WT libraries. The comparisons cover 79 targets, four counting strategies and two editors (WT and HiFi), each compared with NoEP. [Run instructions and workbook verification](manuscript-reads.md#reproduce-the-submitted-calls-from-the-original-fastqs)

## Why the residual at `1450_OT_0235_ALT` leaves the calls unchanged

A positive call requires both Δindel rate (editor minus NoEP) **above 0.1 percentage points** and a BH-adjusted **q below 0.05**. The largest Δ and q differences in the panel occurred at this target under different strategies. Values in the table are rounded.

| Comparison | Submitted Δ (pp) | Fresh Δ (pp) | Submitted q | Fresh q | Call in both analyses |
| --- | ---: | ---: | ---: | ---: | --- |
| WT, Unfiltered | −0.102572 | −0.103837 | 1 | 1 | Negative |
| HiFi, Dedup | 0.074083 | 0.073974 | 0.028681 | 0.028831 | Negative |

For WT Unfiltered, both estimates place the WT indel rate below NoEP. For HiFi Dedup, both q values meet the statistical cutoff, but both Δ estimates are below 0.1 percentage points. All eight editor–strategy comparisons at this target were negative in both analyses.

Across the full panel, the largest absolute differences were 0.001265 percentage points for Δ and 0.000150 for q. The corresponding maxima for the NoEP and edited rates were 0.000133 and 0.001230 percentage points, respectively.

### pBC completeness across all 711 WT batches

The audit compared pBC identities in the archived guide-window and full-length allele tables for each donor, replicate and target (3 donors × 3 replicates × 79 targets). Luminary job `21338827` independently reconstructed the family lists for all 711 batches from the accepted FASTQ run’s target-split read pairs using CRISPResso2 2.3.2 and fastp 0.24.1.

| Evidence available for the membership check | Batches | Result |
| --- | ---: | --- |
| Complete historical CRISPRessoBatch success list | 223 | Both tables matched the success list in 222 batches. At D3908/T1, `1450_OT_0235_ALT`, both omitted the same 932 successful pBCs. |
| FASTQ reconstruction only (complete historical success list unavailable) | 488 | Both tables matched the reconstructed pBC list in all 488 batches. |

Across the other **710 batches**, the audit found no additional missing pBCs. The archived and reconstructed numerators and denominators agreed exactly for both Dedud and Dedud+Dedup in every one of those batches.

### Full-batch pBC check for WT D3908/T1 at `1450_OT_0235_ALT`

The individual CRISPResso outputs matched the replay for all 10,929 successful families. Luminary job `21338660` ran the membership check from commit `fdc16b2` against the two archived tables, comparing the unique pBCs in each table’s `UMI` column with the CRISPRessoBatch success list. The guide-window and full-length tables each contained the same 9,997 successful families and zero unexpected pBCs. The check exited with code 1 as expected and listed **the same 932 missing pBCs** for both tables. Every missing ID matched the earlier family audit.

The 932 families accounted for 17,308 aligned reads, including four indel reads. Eight of their pBCs passed the historical Dedud rule. Their individual allele archives confirmed 2,946 aligned reads and one indel read, which belonged to `AAAGGATAAGGAATTG` (508 aligned reads). All eight representative alleles were unmodified.

| Counting strategy | Archived count | Count with the eight families restored |
| --- | ---: | ---: |
| Dedud (indel reads / aligned reads) | 3 / 20,527 | 4 / 23,473 |
| Dedud+Dedup (indel-positive families / retained families) | 0 / 93 | 0 / 101 |

### Historical family collection

The omissions clustered at the start of the batch. Times below are family completion times recorded on the project disk on April 25, 2025 (UTC+8).

| Family completion time | Omitted from both tables | Included in both tables |
| --- | ---: | ---: |
| 14:56:40–14:59:58 | 929 | 0 |
| 14:59:59–15:00:00 | 3 | 6 |
| 15:00:01–15:33:13 | 0 | 9,991 |

The collector used four batch-level files to mark a target complete. It searched shared memory for each family's guide-window file, then used each match to read both the guide-window table and the full-length allele ZIP. It did not reconcile the matches with the successful pBC list. A missing guide-window file could therefore exclude a family from both tables even when its alignment had succeeded on disk.

### Three-family restart control on Luminary

Job `21338617` used CRISPResso2 2.3.2, fastp 0.24.1 and the historical collector on three WT D3908/T1 families at `1450_OT_0235_ALT`. After a restart with all files present, we moved only the guide-window file for pBC `AAAGGATAAGGAATTG`, leaving its success record and full-length allele ZIP intact.

| Condition | Families reported successful by CRISPRessoBatch | Families in guide-window table | Families in full-length table | Reads in each table |
| --- | ---: | ---: | ---: | ---: |
| All files present, restart with `--no_rerun` | 3/3 | 3 | 3 | 882 |
| Window file removed, same restart command | 3/3 | 2 | 2 | 374 |
| Window file restored, collection only | No new batch run | 3 | 3 | 882 |

With the file missing, `--no_rerun` exited with code 0 and left the window file missing. Both tables omitted the family's 508 reads, including one indel read. Restoring that file alone and running the collector again recovered all 508 reads.

Job `21338648` repeated this control with the [agent’s historical per-pBC collector](manuscript-reads.md#collect-a-historical-per-pbc-batch). With the window file removed, CRISPRessoBatch still reported three successes. The collector exited 1 before counting, named `AAAGGATAAGGAATTG` as missing from both the guide-window and full-length tables, and published neither table. After the file was restored, collection alone exited 0. Both tables contained all three families and 882 reads each. The verified exports matched the historical collector’s output byte for byte at baseline, after an intact restart and after restoration.

**Cause of the 2025 omissions.** On April 27, Slurm job 5010 was cancelled after 74 seconds, and job 5011 used the same script name. The preserved script enabled `--no_rerun`. The control confirms that a missing window file is sufficient for a successful restart to yield incomplete tables. The shared-memory file list from the original collection is unavailable, so the specific operation responsible for the 932 historical omissions remains unknown.

## pBC completeness across all NoEP and HiFi batches

For all 711 NoEP batches and 632 HiFi batches, the accessible archives contained no complete historical list of successful pBC families. In every batch, the pBCs reconstructed from FASTQs matched those in both archived allele tables exactly.

Luminary job `21341762` used the target-split paired FASTQs from nine NoEP libraries and eight HiFi libraries. It checked all 79 targets in each library against both the guide-window and full-length allele tables.

All **2,686 batch-by-strategy comparisons** matched the historical filtered tables in both numerator and denominator. For Dedud, the audit compared indel reads and total retained reads. For Dedud+Dedup, it compared families with an indel-positive representative and total retained families.

## Rerun the full-panel pBC family audit

From the agent checkout, run the command below where the historical archive and the saved FASTQ family replays are readable. It uses Python and pandas and can run locally or in a Slurm allocation. Set the three input paths to your copies and choose a new output directory.

```bash
python validate_family_panel.py \
  --archive-root "$HISTORICAL_PROJECT" \
  --replay-root "$WT_REPLAY" --replay-root "$NOEP_HIFI_REPLAY" \
  --expect-study --output family-audit
```

The historical project supplies `metadata/NovaSeq3_sample_info_new.csv`, the batch plans in `singleUMI_crispresso_input`, the three historical CRISPResso run directories and the original allele tables. The command reads the window and full tables under `singleUMI_crispresso_WTCas9/alleles_all` and `alleles_all_HiFi_NoEP/alleles_all`. It reads the archived Dedud and Dedud+Dedup tables under `target_umi_barcode_table_WTCas9_HiFi_NoEP/{dedudS_filtered_v3,dedudS_filtered_v3_dedup}`.

Point each replay root at the directory containing the `ARM_pool_DONOR_REPLICATE` folders from the FASTQ family reconstruction. Each folder has a `TARGET.tsv.gz` per-family export and its `TARGET.json` receipt. The audit verifies the batch identity, historical read profile, historical Maxi reference digest, sorted pBC membership digest and count totals before comparison. One combined replay root also works. These family-level inputs come from the historical acceptance runs and are needed in addition to the public sequencing data. Use `--workers N` for parallel batch readers within the CPUs allocated to the command.

| Output | What to read |
| --- | --- |
| `batch_summary.tsv` | One row per arm, donor, replicate and target. See missing and unexpected pBC counts for each archived table, plus archived and replayed numerators and denominators. Each delta is replay minus archive. |
| `pbc_differences.tsv` | One row per discrepant pBC. The expected, replay, window and full flags identify which table omitted or added it. The four signed count changes show its contribution even when changes cancel at batch level. Dedud counts reads. Dedud+Dedup counts one representative per retained family. |
| `sources.jsonl` | One record per batch with the selected evidence, historical records found and their status, the expected family-set digest, the replay receipt and the four archive table paths. |
| `summary.json`, `SHA256SUMS` | Overall acceptance, evidence totals by arm, portable result signatures and checksums of the four output files. |

A historical success list qualifies when successes and failures partition the original batch plan and the quantification table confirms the successes. Multiple complete lists must agree. Where a complete list is unavailable, the expected pBCs come from the FASTQ replay. The evidence column records which comparison supports each conclusion.

With `--expect-study`, exit 0 means the expected family sets, full batch summary and per-pBC differences match the verified study fingerprints, including the 932 omissions described above. Without that flag, exit 0 requires zero discrepancies. Exit 1 delivers a complete audit that differs from the selected expectation. Exit 2 indicates an input or execution error. To verify a delivered directory, run `sha256sum -c SHA256SUMS` inside it.

## Export checks

The fresh run produced 2,054 library–target records. Its 8,216-row count table and 632-row model table matched the previously accepted historical analysis byte for byte. The power and figure checks also passed for 1,264 independently simulated power estimates and 158 Figure 3b source rows.

| Fresh export | SHA-256 |
| --- | --- |
| `library_counts.tsv` | `285a86b3b72dc352589ecf4c1df101c3bf885b3080ba14e52ed4168216714018` |
| `panel_comparison.tsv` | `13a906d81fa96b325e9bdc663b8ddc3ac15c373b5861c4c9148d48de5ceaa53f` |
