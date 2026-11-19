"""Create a D3734/T1 manifest from laboratory FASTQs or the public ENA downloads."""
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


def make_manifest(root, output, ena=False):
    with Path(__file__).with_name("d3734_t1_ena.tsv").open() as handle:
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
            candidates = sorted((root / condition).glob(f"{library}_CKDL*_1.fq.gz"))
            if len(candidates) != 1:
                raise ValueError(f"Expected one {library} R1 file in {root / condition}, found {len(candidates)}")
            r1 = candidates[0].absolute()
            stem = r1.name.removesuffix("_1.fq.gz")
            r2 = r1.with_name(f"{stem}_2.fq.gz")
            pair_key = f"{condition}/{stem}_{{mate}}.fq.gz"
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
    parser.add_argument("fastq_root", type=Path, help="Directory containing NoEP/ and HiFi/, or the four ENA FASTQs")
    parser.add_argument("output", type=Path, help="New manifest TSV")
    parser.add_argument("--ena", action="store_true", help="Verify the four ENA files against their archive MD5 values")
    args = parser.parse_args()
    make_manifest(args.fastq_root, args.output, args.ena)
