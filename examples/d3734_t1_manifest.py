"""Create a D3734/T1 or full study panel manifest from laboratory FASTQs or the public ENA downloads."""
import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path, expected_md5=None):
    digest = hashlib.sha256()
    archive_digest = hashlib.md5() if expected_md5 else None
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
            if archive_digest is not None:
                archive_digest.update(block)
    if archive_digest is not None and archive_digest.hexdigest() != expected_md5:
        raise ValueError(f"ENA MD5 mismatch for {path}")
    return digest.hexdigest()


def make_manifest(root, output, ena=False, panel=False):
    catalog = "study_panel_ena.tsv" if panel else "d3734_t1_ena.tsv"
    with Path(__file__).with_name(catalog).open() as handle:
        libraries = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for record in libraries:
        sample = record["sample_name"]
        condition = sample.split("_")[0]
        if ena:
            r1 = (root / Path(record["r1_url"]).name).absolute()
            r2 = (root / Path(record["r2_url"]).name).absolute()
            pair_key = record["run_accession"] + "_{mate}.fastq.gz"
        else:
            library = record["library_id"]
            patterns = [(f"{library}_CKDL*_1.fq.gz", "_{mate}.fq.gz")]
            if condition == "WT":
                donor_replicate = sample.removeprefix("WT_pool_")
                patterns.append((f"*_LVpool_wtCas9_{donor_replicate}_R1.fastq.gz", "_R{mate}.fastq.gz"))
            candidates = [(p, suffix) for pattern, suffix in patterns for p in (root / condition).glob(pattern)]
            if len(candidates) != 1:
                raise ValueError(f"Expected one {library} R1 file in {root / condition}, found {len(candidates)}")
            path, suffix = candidates[0]
            r1 = path.absolute()
            stem = r1.name.removesuffix(suffix.format(mate=1))
            r2 = r1.with_name(stem + suffix.format(mate=2))
            pair_key = f"{condition}/{stem}{suffix}"
        row = dict(dataset=condition, pair_key=pair_key, r1=r1, r2=r2,
                   sha256_r1=sha256(r1, record["md5_r1"] if ena else None),
                   sha256_r2=sha256(r2, record["md5_r2"] if ena else None),
                   sra_row="", sample_name=sample, library_ID=sample, status="matched")
        if ena:
            row["run_accession"] = record["run_accession"]
        rows.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq_root", type=Path, help="Directory containing library folders by condition, or the ENA FASTQs")
    parser.add_argument("output", type=Path, help="New manifest TSV")
    parser.add_argument("--ena", action="store_true", help="Verify the ENA files against their archive MD5 values")
    parser.add_argument("--panel", action="store_true", help="Include all 26 study libraries across NoEP, WT and HiFi")
    args = parser.parse_args()
    make_manifest(args.fastq_root, args.output, args.ena, args.panel)
