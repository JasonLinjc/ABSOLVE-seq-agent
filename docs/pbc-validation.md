# Published pBC annotation validation across 26 libraries

On Luminary, the default published pBC annotation and the annotation rebuilt from the original six Maxi plasmid pool FASTQs with `--published-maxi` produced byte-identical count and model tables **within each read mode**. The historical two-tBC exception remains a separate option (`--historical-maxi`). [Build commands and rule definitions](historical-maxi.md#build-the-reference-once)

| Read mode | Alignments shared by the two annotation runs | Library counts | Model contrasts |
| --- | --- | --- | --- |
| Fast | Archived fast-mode CRISPResso2 alignments | 8,216/8,216 identical | 632/632 identical under each of two donor-weighting schemes |
| Manuscript reproduction | One paired-end CRISPResso2 2.2.13 replay with FLASH | 8,216/8,216 identical | 632/632 identical with equal donor weights |

The 26 libraries comprised nine NoEP, nine WT and eight HiFi libraries. Each annotation run independently counted 79 targets under four strategies. The count comparison covered indel numerators, denominators and pBC family counts. Each model compared WT and HiFi with NoEP for every target and strategy. Control and edited indel rates, Δindel rates, one-sided p values, BH-adjusted q values and positive calls all matched exactly between annotations.

## Manuscript archive check

The original manuscript-mode run retained target-split R1/R2 FASTQs and its default-annotation result tables. Its temporary per-read alignments had been cleared. The replay used the archived FASTQs, the original CRISPResso2 image and the same FLASH parameters. Both recounts read the same 2,054 library/target alignment files and verified their SHA-256 digests. Each model fit weighted the three donors equally on the logit scale.

Both replayed exports (`library_counts.tsv` and `panel_comparison.tsv`) also matched **the original manuscript-mode default archive byte for byte**, including the counts before cross-target filtering. Agreement with the submitted supplementary tables is documented separately. [Supplementary-table comparison](manuscript-reads.md#luminary-validation)

## Export checksums

Each digest below applies to the default and rebuilt annotation exports for that row. The manuscript-mode digests also apply to the archived default exports.

| Read mode | Export | SHA-256 |
| --- | --- | --- |
| Fast | Library counts | `3db1f8203892c29edd231e2d88d942223b9c6a5b44a6d6d23aabc107f90a7a84` |
| Fast | Model contrasts, equal donor weights | `627df9a661d2160184753c39a6bcf39dcef40449c42df7f867b229b1f0342755` |
| Fast | Model contrasts, pBC family weights | `ff94527fbdf0a95bcfc59da57b708ed838692368c34ca69c046c01c537bfaa99` |
| Manuscript | Library counts | `aefad456069a4eb3d8eef4e2cdef9b913c6e8d2791df2be04ea306f60b9fef40` |
| Manuscript | Model contrasts, equal donor weights | `3f52e3068390b9bf5294a816180d2dde7107a9cd29458f0f50b0079341c6290d` |
