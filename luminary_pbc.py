"""Submit or resume the full pBC annotation check and retrieve its verified reports."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from shlex import quote
import subprocess
import sys
import tempfile
import time
from uuid import uuid4


REPORTS = ("summary.tsv", "differences.tsv")
ACTIVE = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED",
          "REQUEUE_HOLD", "RESIZING", "SUSPENDED", "STOPPED", "SIGNALING", "STAGE_OUT"}


def remote(host, command):
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9._@-]*", host):
        raise ValueError("Expected an SSH host alias or user@hostname")
    ssh = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
           "-o", "ControlMaster=no", "-o", "ControlPath=none", host,
           "bash -lc " + quote(command)]
    return subprocess.run(ssh, check=True, stdout=subprocess.PIPE).stdout


def submit(host, repo, dest):
    if not repo:
        repo = remote(host, 'cat "$HOME/.config/absolve-seq/pbc-validation-repo"').decode().strip()
    if not PurePosixPath(repo).is_absolute() or len(repo.splitlines()) != 1:
        raise ValueError("Expected an absolute Luminary checkout path (see --repo)")
    run = str(PurePosixPath(repo, "runs", "pbc-validation-" + uuid4().hex))
    command = (f"cd -- {quote(repo)} && test -r data/pbc-validation.json && "
               f"sbatch --parsable -p cpu1,cpu2 --output=slurm-%j.out "
               f"validate_pbc_panel.sbatch --output {quote(run)}")
    reply = remote(host, command).decode().strip()
    job = reply.partition(";")[0]
    if not job.isdigit():
        raise ValueError(f"Unexpected sbatch response: {reply!r}. Check squeue before submitting again.")
    receipt = dict(host=host, job=job, remote_run=run, log=f"{repo}/slurm-{job}.out")
    (dest / "job.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def wait_for_job(host, job):
    if not re.fullmatch(r"[0-9]+", job):
        raise ValueError("Invalid Slurm job ID")
    command = f"sacct -X -n -P -j {job} -o JobIDRaw,State,ExitCode"
    previous = None
    while True:
        rows = remote(host, command).decode().splitlines()
        for row in rows:
            job_id, state, code = row.split("|")
            if job_id == job:
                state = state.split()[0]  # CANCELLED can include the cancelling user.
                if state != previous:
                    print(f"Slurm {job} {state} ({code})", flush=True)
                    previous = state
                if state not in ACTIVE:
                    return state, code
        time.sleep(20)


def collect(receipt, dest, state, code):
    host, run = receipt["host"], receipt["remote_run"]
    prefix = f"cd -- {quote(run)} && "
    expected = {("COMPLETED", "0:0"): "pass", ("FAILED", "1:0"): "different"}.get((state, code))
    if not expected:
        raise ValueError(f"Job ended with {state} ({code}). See {receipt['log']}")
    status = json.loads(remote(host, prefix + "cat run.json"))
    if status["status"] != expected:
        raise ValueError(f"Unexpected validation status {status['status']!r}. See {receipt['log']}")
    manifest = remote(host, prefix + "sha256sum --check --strict --status SHA256SUMS && cat SHA256SUMS")
    digests = {}
    for line in manifest.decode().splitlines():
        digest, name = line.split("  ", 1)
        digests[name] = digest
    with tempfile.TemporaryDirectory(dir=dest) as tmp:
        stage = Path(tmp)
        for name in REPORTS:
            data = remote(host, prefix + "cat " + name)
            if hashlib.sha256(data).hexdigest() != digests[name]:
                raise ValueError(f"Downloaded digest mismatch for {name}")
            (stage / name).write_bytes(data)
        (stage / "SHA256SUMS").write_bytes(manifest)
        for name in (*REPORTS, "SHA256SUMS"):
            (stage / name).replace(dest / name)
    print(f"{expected.upper()} (differing fields or file checks: {status['differences']}). Reports in {dest}")
    return int(expected == "different")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--host", default="luminary", help="SSH alias (default: luminary)")
    cli.add_argument("--repo", help="Absolute configured repository path on Luminary")
    group = cli.add_mutually_exclusive_group()
    group.add_argument("--output", type=Path, help="New local delivery directory")
    group.add_argument("--resume", type=Path, help="Local directory containing a previous job.json")
    args = cli.parse_args()
    if args.resume:
        dest = args.resume.resolve()
        receipt = json.loads((dest / "job.json").read_text())
    else:
        dest = (args.output or Path("runs", "pbc-delivery-" + uuid4().hex)).resolve()
        dest.mkdir(parents=True, exist_ok=False)
        receipt = submit(args.host, args.repo, dest)
    print(f"Job {receipt['job']} on {receipt['host']}. Resume with --resume {quote(str(dest))}", flush=True)
    state, code = wait_for_job(receipt["host"], receipt["job"])
    return collect(receipt, dest, state, code)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        print(f"Validation delivery error: {exc}", file=sys.stderr)
        sys.exit(2)
