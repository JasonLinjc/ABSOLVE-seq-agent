"""Create a D3734/T1 manifest from the original, compressed NoEP and HiFi FASTQs."""
import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def make_manifest(root, output):
    rows = []
    for condition, library in (("NoEP", "LV_S3"), ("HiFi", "LV_S1")):
        candidates = sorted((root / condition).glob(f"{library}_CKDL*_1.fq.gz"))
        if len(candidates) != 1:
            raise ValueError(f"Expected one {library} R1 file in {root / condition}, found {len(candidates)}")
        r1 = candidates[0].absolute()
        stem = r1.name.removesuffix("_1.fq.gz")
        r2 = r1.with_name(f"{stem}_2.fq.gz")
        sample = f"{condition}_pool_D3734_T1"
        rows.append(dict(dataset=condition, pair_key=f"{condition}/{stem}_{{mate}}.fq.gz",
                         r1=r1, r2=r2, sha256_r1=sha256(r1), sha256_r2=sha256(r2),
                         sra_row="", sample_name=sample, library_ID=sample, status="matched"))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq_root", type=Path, help="Directory containing NoEP/ and HiFi/")
    parser.add_argument("output", type=Path, help="New manifest TSV")
    args = parser.parse_args()
    make_manifest(args.fastq_root, args.output)
