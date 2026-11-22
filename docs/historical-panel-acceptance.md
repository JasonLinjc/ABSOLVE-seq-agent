# Acceptance of the historical panel from raw FASTQs

The Luminary run reproduced **all 632 classifications in Supplementary Tables 2–5 (63 positive and 569 negative)** from 52 checksum-verified raw FASTQs. It processed 334,970,765 read pairs from 26 libraries, including all nine WT libraries. The comparisons cover 79 targets, four counting strategies and two editors (WT and HiFi), each compared with NoEP. [Run instructions and workbook verification](manuscript-reads.md#reproduce-the-submitted-calls-from-the-original-fastqs)

## Why the residual at `1450_OT_0235_ALT` leaves the calls unchanged

A positive call requires both Δindel rate (editor minus NoEP) **above 0.1 percentage points** and a BH-adjusted **q below 0.05**. The largest Δ and q differences in the panel occurred at this target under different strategies. Values in the table are rounded.

| Comparison | Submitted Δ (pp) | Fresh Δ (pp) | Submitted q | Fresh q | Call in both analyses |
| --- | ---: | ---: | ---: | ---: | --- |
| WT, Unfiltered | −0.102572 | −0.103837 | 1 | 1 | Negative |
| HiFi, Dedup | 0.074083 | 0.073974 | 0.028681 | 0.028831 | Negative |

For WT Unfiltered, both estimates place the WT indel rate below NoEP. For HiFi Dedup, both q values meet the statistical cutoff, but both Δ estimates are below 0.1 percentage points. All eight editor–strategy comparisons at this target were negative in both analyses.

The family audit traced a count difference to the historical combined allele tables for WT donor D3908, replicate T1. Eight pBC families present in the individual alignment outputs and retained by Dedud were absent from the combined tables. Including them added 2,946 reads, one carrying an indel. Across the full panel, the largest absolute differences were 0.001265 percentage points for Δ and 0.000150 for q. The corresponding maxima for the NoEP and edited rates were 0.000133 and 0.001230 percentage points, respectively.

## Export checks

The fresh run produced 2,054 library–target records. Its 8,216-row count table and 632-row model table matched the previously accepted historical analysis byte for byte. The power and figure checks also passed for 1,264 independently simulated power estimates and 158 Figure 3b source rows.

| Fresh export | SHA-256 |
| --- | --- |
| `library_counts.tsv` | `285a86b3b72dc352589ecf4c1df101c3bf885b3080ba14e52ed4168216714018` |
| `panel_comparison.tsv` | `13a906d81fa96b325e9bdc663b8ddc3ac15c373b5861c4c9148d48de5ceaa53f` |
