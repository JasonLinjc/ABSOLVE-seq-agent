# Maxi pBC reference builds and manuscript results

Use the published pBC annotation directly for the default analysis. To rebuild qualifying pBCs from the original study’s six Maxi FASTQs, choose a rule. Both rules pool observations across all three libraries, require the designed insert to match in every observation for a given pBC/tBC, and exclude pBCs associated with more than two observed tBCs.

| Rule | One observed tBC | Two observed tBCs | Builder | Analysis method |
| --- | --- | --- | --- | --- |
| Published | Entropy > 1.477 bits | Entropy > 1.477 bits | `--published-maxi` | `manuscript` |
| Historical | Entropy > 1.477 bits | Any entropy | `--historical-maxi` | `historical` |

## Build the reference once

Use the original six compressed Maxi FASTQs listed in [the input manifest](../runtime/maxi/inputs.sha256). The hashes identify the original compressed files. An SRA export with different FASTQ headers or gzip encoding can have a different hash. For the build, use Java with a Trimmomatic 0.36 jar, Python 3.8 or later, Bash, awk and GNU coreutils on a Linux compute node. The working filesystem needs space for the trimmed R1s and sort files. On Luminary, allow at least 130 GiB of free scratch and 32 GiB of RAM for the following configuration.

```bash
# Rebuild the published rule.
python setup_references.py --published-maxi /data/original_maxi \
  --trimmomatic-jar /software/trimmomatic-0.36.jar --threads 15 \
  --output /data/references/published-maxi

# Alternatively, rebuild the historical rule.
python setup_references.py --historical-maxi /data/original_maxi \
  --trimmomatic-jar /software/trimmomatic-0.36.jar --threads 15 \
  --output /data/references/historical-maxi
```

The builder verifies all six input hashes and each of the 79 target memberships against the selected rule’s SHA-256 manifest before publishing the new directory. The published-rule manifest was computed independently from the pinned published annotation by taking the first row for each target/pBC, selecting qualifying barcodes, and hashing the sorted barcode list for each target. Only counts and digests from that comparison are packaged with the builder. Outputs include one sorted pBC list per target, a combined `target_pbc.tsv.gz`, per-target counts and hashes in `target_comparison.tsv`, trimming logs and `provenance.json`. The published rule has **5,190,216 target/pBC keys** and the historical rule has **5,248,941**. Reference files are checked again when read by a worker. Keep the directory accessible from all worker nodes and bind it into the container when using a path outside the standard mounts. The [legacy fast guide](manuscript-reads.md#legacy-fast-compatibility) describes the older pool-voting methods.

## Analyze a matched pair or the panel

Prepare the normal sample manifest and target definitions as in the [quick start](d3734-t1-quick-start.md). Configure manuscript read-processing dependencies as in [the read-mode guide](manuscript-reads.md). Use an absolute reference path and a new run directory. To reuse a build, select the method recorded in its `provenance.json`. Only the first build needs the Maxi FASTQs.

```bash
# Published-rule rebuild, using the default unified paired-end workflow.
python paired.py --donor D3734 --replicate T1 \
  --maxi-method manuscript --maxi-reference /data/references/published-maxi \
  --output runs/D3734-T1-rebuilt-published

# Historical-rule examples using manuscript read processing.
python paired.py --read-mode manuscript --donor D3734 --replicate T1 \
  --maxi-method historical --maxi-reference /data/references/historical-maxi \
  --output runs/D3734-T1-historical

python paired.py --manuscript-panel --read-mode manuscript \
  --maxi-method historical --maxi-reference /data/references/historical-maxi \
  --output runs/panel-historical
```

The paired output has 316 target/strategy records. The full model estimates 632 editor/target/strategy contrasts. The new reference governs Dedud and Dedud+Dedup. The read mode and its allele tie handling are separate choices from the Maxi rule. For existing fast-mode comparisons, add `--read-mode fast` when submitting the pair with `sbatch paired.sbatch`.

## Refit existing counts

A sufficient-count refit preserves the read processing, allele tie handling and pBC filtering already applied to the input counts. To change any of those choices, run the corresponding read or family analysis first. The refit checks all 8,216 library/target/strategy keys, fits the equal-donor model and records the count-table hash. Use the same reference-data setup as the manuscript panel. Provide R with brglm2 or set `ABSOLVE_IMAGE` to the statistics image.

```bash
python paired.py --manuscript-panel --refit-counts previous/results/library_counts.tsv \
  --output runs/panel-refit
```

## Generate Figure 3b

Follow the [Figure 3b guide](manuscript-figures.md) to export the Dedud plot and its Excel source data.

## Luminary acceptance

Using the original six Maxi FASTQs with `--published-maxi` recovered all **5,190,216** published target/pBC keys. An independent set comparison against the published annotation found zero extra or missing keys for each of the 79 targets. A D3734/T1 fast pair using only the rebuilt pBC reference produced the same 316-row `paired_target_rates.tsv` as the default annotation, byte for byte (SHA-256 `0b0465ecc5aa244cc2a10654e153dbe8467322cbedd4f5bfc72806058018980b`). Both runs reused the same checksum-verified NoEP and HiFi alignments and independently counted all four strategies. The historical-rule pair also reproduced its earlier result.

### Published-rule rebuild across 26 libraries

[Full-panel validation in fast and manuscript reproduction modes, including the manuscript archive check and export checksums](pbc-validation.md)

### Historical-rule rebuild and saved-count refit

The historical-rule build read the original six FASTQs and recovered 5,248,941 keys. All 79 per-target SHA-256 values matched the historical archive. The combined export reproduced the independently measured per-target comparison table. A fresh D3734/T1 NoEP–HiFi manuscript-mode run used this rebuilt membership. All 316 exported rows matched independent counts from the archived-label replay, including both arm numerators, denominators, rates and every defined difference.

A separate 26-library model refit used the existing encounter-order analysis counts. The model-only command validated all 8,216 count cells and reproduced the preceding 632-row results file byte for byte. Compared with the submitted supplements, all 632 primary calls agree. The largest absolute Δ and q differences were 0.001265 percentage points and 0.000150, respectively. The two-library FASTQ test and the 26-library count refit exercise different stages.
