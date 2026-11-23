"""Audit historical pBC completeness across the 26 study libraries and 79 targets."""
import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory

from panel_report.family_audit import METRICS, PRESENCE, RUN_DIRS, audit_batch

TARGETS = Path(__file__).resolve().parent / "runtime/maxi/reference_sha256.tsv"
# Portable fingerprints of the independently cross-checked historical audit.
STUDY = {
    "batch_summary.tsv": "19c936c12db4367db787b67033b5843f9ffed0ddd8dea480795ef48ec16c41f1",
    "pbc_differences.tsv": "74b9278c47302594eb197d870831865fdbd40a86fedfe164a0444af177335474",
    "expected_family_sets": "9537ef347f2bf5c9536d5e54e740b0b2886efc7410e15c7a0b531ddfacb15d01",
}
ARMS = {"noEP": "NoEP", "HiFiCas9": "HiFi", "WTCas9": "WT"}


def panel_tasks(archive, roots):
    for directory in (*RUN_DIRS, "singleUMI_crispresso_input"):
        if not (archive / directory).is_dir():
            raise ValueError(f"Missing historical record search directory: {archive / directory}")
    targets = [r["target"] for r in csv.DictReader(TARGETS.read_text().splitlines(), delimiter="\t")]
    if len(set(targets)) != 79 or len(targets) != 79:
        raise ValueError(f"Expected 79 unique targets in {TARGETS}")
    ref_hash = hashlib.sha256(TARGETS.read_bytes()).hexdigest()
    with (archive / "metadata/NovaSeq3_sample_info_new.csv").open() as handle:
        libraries = sorted((ARMS[r["cas_type"]], r["donor_id"], r["replicate"], r["LV_ID"])
                           for r in csv.DictReader(handle) if r["cas_type"] in ARMS)
    expected = {(a, d, r) for a in ARMS.values() for d in ("D3734", "D3742", "D3908")
                for r in ("T1", "T2", "T3")} - {("HiFi", "D3734", "T3")}
    if {lib[:3] for lib in libraries} != expected or len({lib[3] for lib in libraries}) != 26 or len(libraries) != 26:
        raise ValueError("Sample metadata do not identify the 26 distinct historical study libraries")
    for arm, donor, rep, lv in libraries:
        name = f"{arm}_pool_{donor}_{rep}"
        matches = [root / name for root in roots if (root / name).is_dir()]
        if len(matches) != 1:
            raise ValueError(f"Expected one replay directory for {name}, found {matches}")
        for target in sorted(targets):
            yield archive, matches[0], arm, donor, rep, target, lv, ref_hash


def write_audit(staged, tasks, workers, expect_study=False):
    family_sets = hashlib.sha256()
    counts = Counter()
    evidence = Counter()
    details = ["arm", "donor", "replicate", "target", "evidence", "pbc", *PRESENCE]
    details += [label + "_" + metric for metric in METRICS for label in ("archived", "replay", "delta")]
    with ExitStack() as stack:
        def output(name):
            return stack.enter_context((staged / name).open("w", newline=""))
        batch = output("batch_summary.tsv")
        sources = output("sources.jsonl")
        detail = csv.DictWriter(output("pbc_differences.tsv"), fieldnames=details, delimiter="\t", lineterminator="\n")
        detail.writeheader()
        executor = stack.enter_context(ProcessPoolExecutor(workers)) if workers > 1 else None
        results = executor.map(audit_batch, tasks, chunksize=4) if executor else map(audit_batch, tasks)
        batches = None
        for row, differences, source in results:
            if batches is None:
                batches = csv.DictWriter(batch, fieldnames=list(row), delimiter="\t", lineterminator="\n")
                batches.writeheader()
            batches.writerow(row)
            detail.writerows(differences)
            sources.write(json.dumps(source, sort_keys=True) + "\n")
            fields = ("arm", "donor", "replicate", "target", "expected_ids_sha256")
            family_sets.update(("\t".join(source[k] for k in fields) + "\n").encode())
            evidence[(row["arm"], row["evidence"])] += 1
            counts.update(batches=1, discrepant_batches=int(bool(differences)), discrepant_pbcs=len(differences))
    summary = dict(counts, evidence=[dict(arm=a, evidence=e, batches=n) for (a, e), n in sorted(evidence.items())])
    signatures = {name: hashlib.sha256((staged / name).read_bytes()).hexdigest()
                  for name in ("batch_summary.tsv", "pbc_differences.tsv")}
    signatures["expected_family_sets"] = family_sets.hexdigest()
    summary["signatures"] = signatures
    if expect_study:
        summary["baseline_mismatches"] = [k for k in STUDY if signatures[k] != STUDY[k]]
        summary["accepted"] = not summary["baseline_mismatches"]
    else:
        summary["accepted"] = counts["discrepant_batches"] == 0
    (staged / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    files = sorted(staged.iterdir())
    with (staged / "SHA256SUMS").open("w") as out:
        for path in files:
            out.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", required=True, type=Path, help="Root of the original historical project")
    parser.add_argument("--replay-root", required=True, action="append", type=Path,
                        help="Directory of receipt-checked FASTQ family replays; repeat for separate roots")
    parser.add_argument("--output", type=Path, required=True, help="New directory for this acceptance run")
    parser.add_argument("--workers", type=int, default=1, help="Parallel batch readers (default 1)")
    parser.add_argument("--expect-study", action="store_true",
                        help="Accept the exact study audit, including its documented historical omission")
    args = parser.parse_args()
    try:
        if args.workers < 1 or args.output.exists():
            raise ValueError("Use a positive worker count and a new output directory")
        tasks = list(panel_tasks(args.archive_root.resolve(), [p.resolve() for p in args.replay_root]))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=args.output.parent) as tmp:
            staged = Path(tmp) / "audit"
            staged.mkdir()
            summary = write_audit(staged, tasks, args.workers, args.expect_study)
            staged.rename(args.output)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"Family audit incomplete: {error}\n")
    print(f"{summary['batches']} batches; {summary['discrepant_batches']} with differences. Results: {args.output}")
    print("Acceptance:", "PASS" if summary["accepted"] else "DIFFERENCES")
    return int(not summary["accepted"])


if __name__ == "__main__":
    sys.exit(main())
