# Historical Maxi reference and results figures

The historical rule pools the three original Maxi libraries before assigning a target/pBC key. A pBC/tBC needs exact agreement with its designed insert in every observation. A pBC with one observed tBC needs entropy above 1.477 bits. A pBC with exactly two observed tBCs can pass at any entropy. Count all distinct observed tBCs when applying this rule. The default `manuscript` method uses the upstream published annotation and applies its entropy cutoff to both barcode counts.

## Build the reference once

Use the original six compressed Maxi FASTQs listed in [the input manifest](../runtime/maxi/inputs.sha256). The hashes identify the original compressed files. An SRA export with different FASTQ headers or gzip encoding can have a different hash. For the build, use Java with a Trimmomatic 0.36 jar, Python 3.8 or later, Bash, awk and GNU coreutils on a Linux compute node. The working filesystem needs space for the trimmed R1s and sort files. On Luminary, allow at least 130 GiB of free scratch and 32 GiB of RAM for the following configuration.

```bash
python setup_references.py --historical-maxi /data/original_maxi \
  --trimmomatic-jar /software/trimmomatic-0.36.jar --threads 15 \
  --output /data/references/historical-maxi
```

The builder verifies all six input hashes and each of the 79 target memberships against the archived SHA-256 manifest before publishing the new directory. Outputs include one sorted pBC list per target, a combined `target_pbc.tsv.gz`, per-target counts and hashes in `target_comparison.tsv`, trimming logs and `provenance.json`. The expected total is **5,248,941 target/pBC keys**. Reference files are checked again when read by a worker. Keep the directory accessible from all worker nodes and bind it into the container when using a path outside the standard mounts. See the main README for the alternative pool-voting methods.

## Analyze a matched pair or the panel

Prepare the normal sample manifest and target definitions as in the [quick start](d3734-t1-quick-start.md). Configure manuscript read-processing dependencies as in [the read-mode guide](manuscript-reads.md). Use an absolute reference path and a new run directory.

```bash
python paired.py --read-mode manuscript --donor D3734 --replicate T1 \
  --maxi-method historical --maxi-reference /data/references/historical-maxi \
  --output runs/D3734-T1-historical

python paired.py --manuscript-panel --read-mode manuscript \
  --maxi-method historical --maxi-reference /data/references/historical-maxi \
  --output runs/panel-historical
```

The paired output has 316 target/strategy records. The full model estimates 632 editor/target/strategy contrasts. The new reference governs Dedud and Dedud+Dedup. The read mode and its allele tie handling are separate choices from the Maxi rule. A fast-mode pair also accepts the historical reference via `sbatch paired.sbatch`.

## Refit existing counts

A sufficient-count refit preserves the read processing, allele tie handling and pBC filtering already applied to the input counts. To change any of those choices, run the corresponding read or family analysis first. The refit checks all 8,216 library/target/strategy keys, fits the equal-donor model and records the count-table hash. Use the same reference-data setup as the manuscript panel. Provide R with brglm2 or set `ABSOLVE_IMAGE` to the statistics image.

```bash
python paired.py --manuscript-panel --refit-counts previous/results/library_counts.tsv \
  --output runs/panel-refit
```

## Generate figures

The [manuscript figure guide](manuscript-figures.md) covers Figures 2a, 2b, 2d, 3b and Supplementary Figures 10–13. Figure 2d uses the original published strict entropy index. The historical two-tBC exception is a separate analysis rule.

## Luminary acceptance

The production reference builder read the original six FASTQs and recovered 5,248,941 keys. All 79 per-target SHA-256 values matched the historical archive. The combined export reproduced the independently measured per-target comparison table. A fresh D3734/T1 NoEP–HiFi manuscript-mode run used this rebuilt membership. All 316 exported rows matched independent counts from the archived-label replay, including both arm numerators, denominators, rates and every defined difference.

The 26-library acceptance used the existing encounter-order analysis counts. The model-only command validated all 8,216 count cells and reproduced the preceding 632-row results file byte for byte. Compared with the submitted supplements, all 632 primary calls agree. The largest absolute Δ and q differences were 0.001265 percentage points and 0.000150, respectively. The two-library FASTQ test and the 26-library count refit exercise different stages.
