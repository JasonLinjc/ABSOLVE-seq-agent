# Acceptance of the historical panel from raw FASTQs

The Luminary run reproduced **all 632 classifications in Supplementary Tables 2–5 (63 positive and 569 negative)** from 52 checksum-verified raw FASTQs. It processed 334,970,765 read pairs from 26 libraries, including all nine WT libraries. The comparisons cover 79 targets, four counting strategies and two editors (WT and HiFi), each compared with NoEP. [Run instructions and workbook verification](manuscript-reads.md#reproduce-the-submitted-calls-from-the-original-fastqs)

## Why the residual at `1450_OT_0235_ALT` leaves the calls unchanged

A positive call requires both Δindel rate (editor minus NoEP) **above 0.1 percentage points** and a BH-adjusted **q below 0.05**. The largest Δ and q differences in the panel occurred at this target under different strategies. Values in the table are rounded.

| Comparison | Submitted Δ (pp) | Fresh Δ (pp) | Submitted q | Fresh q | Call in both analyses |
| --- | ---: | ---: | ---: | ---: | --- |
| WT, Unfiltered | −0.102572 | −0.103837 | 1 | 1 | Negative |
| HiFi, Dedup | 0.074083 | 0.073974 | 0.028681 | 0.028831 | Negative |

For WT Unfiltered, both estimates place the WT indel rate below NoEP. For HiFi Dedup, both q values meet the statistical cutoff, but both Δ estimates are below 0.1 percentage points. All eight editor–strategy comparisons at this target were negative in both analyses.

For WT donor D3908, replicate T1, the individual CRISPResso outputs matched the replay for all 10,929 aligned pBC families. The historical combined allele tables contained 9,997 of them. The 932 families omitted during aggregation accounted for 17,308 aligned reads, including four indel reads.

Only eight of those 932 families had pBCs that qualified under the historical Dedud rule. Including them added 2,946 reads and one indel read, changing the Dedud count from 3/20,527 to 4/23,473 (indel reads / aligned reads). All eight families had an unmodified representative allele, so the Dedud+Dedup count changed from 0/93 to 0/101 (indel-positive families / retained families).

Across the full panel, the largest absolute differences were 0.001265 percentage points for Δ and 0.000150 for q. The corresponding maxima for the NoEP and edited rates were 0.000133 and 0.001230 percentage points, respectively.

### Historical family collection

The omissions clustered at the start of the batch. Times below are family completion times recorded on the project disk on April 25, 2025 (UTC+8).

| Family completion time | Omitted from both tables | Included in both tables |
| --- | ---: | ---: |
| 14:56:40–14:59:58 | 929 | 0 |
| 14:59:59–15:00:00 | 3 | 6 |
| 15:00:01–15:33:13 | 0 | 9,991 |

The collector used four batch-level files to mark a target complete. It searched shared memory for each family's guide-window file, then used each match to read both the guide-window table and the full-length allele ZIP. It did not reconcile the matches with the successful pBC list. A missing guide-window file could therefore exclude a family from both tables even when its alignment had succeeded on disk.

**Restart hypothesis.** On April 27, Slurm job 5010 was cancelled after 74 seconds, and job 5011 used the same script name. The preserved script enabled `--no_rerun`, which could skip a previously completed family without checking its guide-window file. The shared-memory file list from collection time is unavailable, so the records cannot establish whether copying, restart behavior or cleanup caused the omissions.

## Export checks

The fresh run produced 2,054 library–target records. Its 8,216-row count table and 632-row model table matched the previously accepted historical analysis byte for byte. The power and figure checks also passed for 1,264 independently simulated power estimates and 158 Figure 3b source rows.

| Fresh export | SHA-256 |
| --- | --- |
| `library_counts.tsv` | `285a86b3b72dc352589ecf4c1df101c3bf885b3080ba14e52ed4168216714018` |
| `panel_comparison.tsv` | `13a906d81fa96b325e9bdc663b8ddc3ac15c373b5861c4c9148d48de5ceaa53f` |
