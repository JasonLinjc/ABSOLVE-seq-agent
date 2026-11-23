"""Compare published and rebuilt pBC annotations on the full retained study panel."""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory, mkdtemp

from absolve.maxi_qc import sha256
from panel_report.manuscript_reads import save, table
from panel_report.validation_counts import AGENT, SOURCES, recount

PATHS = ("data", "fast_inventory", "manuscript_alignments", "manuscript_archive",
         "rebuilt_reference", "submitted_tables", "fast_stats_image", "manuscript_stats_image")
COUNT_KEY = ("target", "strategy", "donor", "replicate", "arm")
MODEL_KEY = ("target", "strategy", "editor")
DIMENSIONS = (*COUNT_KEY, "editor")
DIFF_FIELDS = ("mode", "comparison", "export", *DIMENSIONS, "field", "left", "right")
FAST_BENCHMARKS = {
    "counts": "3db1f8203892c29edd231e2d88d942223b9c6a5b44a6d6d23aabc107f90a7a84",
    "equal-donor": "627df9a661d2160184753c39a6bcf39dcef40449c42df7f867b229b1f0342755",
    "family-rows": "ff94527fbdf0a95bcfc59da57b708ed838692368c34ca69c046c01c537bfaa99",
}


def load_paths(config):
    paths = json.loads(config.read_text())
    if set(paths) != set(PATHS):
        raise ValueError(f"Config requires these path keys: {', '.join(PATHS)}")
    return {k: (config.parent / Path(paths[k]).expanduser()).resolve(strict=True) for k in PATHS}


def fit(stage, counts, cfg, mode, margin):
    stage.mkdir()
    os.link(counts, stage/"library_counts.tsv")
    script = AGENT/"panel_report/manuscript_model.R"
    binds = ",".join(map(str, (stage, cfg["submitted_tables"], script.parent)))
    cmd = ["apptainer", "exec", "--cleanenv", "--bind", binds, cfg[mode+"_stats_image"],
           "Rscript", script, stage, cfg["submitted_tables"], margin]
    with (stage/"model.log").open("w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                              env=dict(os.environ, APPTAINERENV_OPENBLAS_NUM_THREADS="1", APPTAINERENV_OMP_NUM_THREADS="1"))
    if proc.returncode:
        raise RuntimeError(f"Model fit failed; see {stage/'model.log'}")


def keyed(path, keys):
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames
        if not fields or not set(keys) <= set(fields) or len(fields) != len(set(fields)):
            raise ValueError(f"Invalid header: {path}")
        rows = {}
        for row in reader:
            key = tuple(row[k] for k in keys)
            if key in rows or None in row or None in row.values():
                raise ValueError(f"Duplicate key or malformed row in {path}: {key}")
            rows[key] = row
        return fields, rows


def compare(mode, name, export, left, right, differences):
    keys, expected = (COUNT_KEY, 8216) if export == "counts" else (MODEL_KEY, 632)
    a_cols, a = keyed(left, keys)
    b_cols, b = keyed(right, keys)
    start = len(differences)

    def record(key, field, x, y):
        differences.append(dict(mode=mode, comparison=name, export=export,
                                **dict(zip(keys, key)), field=field, left=x, right=y))

    if len(a) != expected or len(b) != expected:
        record((), "rows_expected_" + str(expected), len(a), len(b))
    if a_cols != b_cols:
        record((), "header", ",".join(a_cols), ",".join(b_cols))
    for key in sorted(a.keys() | b.keys()):
        if key not in a or key not in b:
            record(key, "row_present", key in a, key in b)
            continue
        for col in a_cols:
            if col in b_cols and a[key][col] != b[key][col]:
                record(key, col, a[key][col], b[key][col])
    hashes = sha256(left), sha256(right)
    if hashes[0] != hashes[1] and len(differences) == start:
        record((), "byte_sha256", *hashes)
    if mode == "fast" and hashes[0] != FAST_BENCHMARKS[export]:
        record((), "accepted_vs_default_sha256", FAST_BENCHMARKS[export], hashes[0])
    return dict(mode=mode, comparison=name, export=export, left_rows=len(a), right_rows=len(b),
                differences=len(differences)-start, left_sha256=hashes[0], right_sha256=hashes[1])


def validate(cfg, out, workers):
    with TemporaryDirectory(dir=out, prefix="fast-workers-") as scratch:
        counts = recount(cfg, Path(scratch), workers)
    summaries, differences = [], []
    exports = {}
    for mode in ("fast", "manuscript"):
        margins = ("equal-donor", "family-rows") if mode == "fast" else ("equal-donor",)
        for source in SOURCES:
            folder = out/mode/source
            folder.mkdir(parents=True)
            count_file = folder/"library_counts.tsv"
            save(count_file, counts[mode, source])
            exports[mode, source, "counts"] = count_file
            for margin in margins:
                stage = folder/margin
                print(f"Fitting {mode}/{source}/{margin}", flush=True)
                fit(stage, count_file, cfg, mode, margin)
                exports[mode, source, margin] = stage/"panel_comparison.tsv"
        for key in ("counts", *margins):
            ref = exports[mode, "default", key]
            summaries.append(compare(mode, "default_vs_rebuilt", key, ref, exports[mode, "rebuilt", key], differences))
            if mode == "manuscript":
                snapshot = out/mode/"archive"/ref.name
                snapshot.parent.mkdir(exist_ok=True)
                shutil.copyfile(cfg["manuscript_archive"]/ref.name, snapshot)
                exports[mode, "archive", key] = snapshot
                summaries.append(compare(mode, "archive_vs_default", key, snapshot, ref, differences))
    save(out/"summary.tsv", summaries)
    with (out/"differences.tsv").open("x", newline="") as f:
        writer = csv.DictWriter(f, DIFF_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(differences)
    with (out/"SHA256SUMS").open("x") as f:
        for path in sorted(set(exports.values()) | {out/"summary.tsv",out/"differences.tsv"}):
            f.write(f"{sha256(path)}  {path.relative_to(out)}\n")
    return len(differences)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--config", type=Path, default=Path("data/pbc-validation.json"))
    cli.add_argument("--output", type=Path, help="New directory (default: a unique directory under runs/)")
    cli.add_argument("--workers", type=int, default=int(os.getenv("SLURM_CPUS_PER_TASK", "1")))
    args = cli.parse_args()
    if args.workers < 1:
        cli.error("--workers must be positive")
    cfg = load_paths(args.config.resolve())
    if args.output:
        out = args.output.resolve()
        out.mkdir(parents=True, exist_ok=False)
    else:
        runs = Path("runs").resolve()
        runs.mkdir(exist_ok=True)
        out = Path(mkdtemp(prefix="pbc-validation-", dir=runs))
    print(f"Results: {out}", flush=True)
    receipt = dict(paths={k:str(v) for k,v in cfg.items()}, workers=args.workers, status="running")
    path = out/"run.json"
    path.write_text(json.dumps(receipt, indent=2)+"\n")
    try:
        mismatches = validate(cfg, out, args.workers)
    except Exception as exc:
        receipt.update(status="error", error=str(exc))
        raise
    else:
        receipt.update(status="different" if mismatches else "pass", differences=mismatches)
    finally:
        path.write_text(json.dumps(receipt, indent=2)+"\n")
    print(f"{receipt['status'].upper()}: {out/'summary.tsv'}", flush=True)
    return int(bool(mismatches))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Validation error: {error}", file=sys.stderr)
        sys.exit(2)
