"""Recount the original panel from retained fast and paired-end alignments."""
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from absolve.historical_maxi import rebuilt_labels, reference, targets
from absolve.maxi_qc import PUBLISHED_SHA256, sha256
from panel_report.manuscript_reads import STRATEGIES, read_aligned_families, summarize, table
from absolve.target_identity import TargetIdentity

AGENT = Path(__file__).resolve().parents[1]
SOURCES = ("default", "rebuilt")
ARMS = {"WTCas9": "WT", "HiFiCas9": "HiFi", "noEP": "NoEP"}


def prepare(cfg, work):
    refs = cfg["data"] / "references"
    if (refs/"published/source.sha256").read_text().strip() != PUBLISHED_SHA256:
        raise ValueError("Published pBC annotation source digest differs")
    rebuilt = reference(cfg["rebuilt_reference"], "manuscript")
    design = {f"{ARMS[r['cas_type']]}_pool_{r['donor_id']}_{r['replicate']}": r
              for r in table(refs/"target_info/NovaSeq3_sample_info.csv", ",")}
    fast = table(cfg["fast_inventory"])
    manuscript = table(cfg["manuscript_archive"]/"sample_inventory.tsv")
    expected = {r["sample_name"]: r for r in manuscript}
    if len(design) != 26 or len(fast) != 26 or len(manuscript) != 26 or not (
            set(design) == set(expected) == {r["sample"] for r in fast}):
        raise ValueError("Expected the same 26 distinct study libraries in both inventories")
    for s in fast:
        name = s["sample"]
        if name != f"{s['arm']}_pool_{s['donor']}_{s['replicate']}" or any(
                s[k] != expected[name][k] for k in ("donor", "replicate", "arm", "sha256_r1", "sha256_r2")):
            raise ValueError(f"Library identity differs between inventories: {name}")
        folder = Path(s["folder"]).resolve(strict=True)
        prep = json.loads((folder/"prepared/summary.json").read_text())
        sample_id = name if s["arm"] == "WT" else design[name]["LV_ID"]
        if not (folder.name == s["arm"] and prep["sample_id"] == sample_id and
                prep["scope"] == "full" and prep["status"] == "complete" and prep["counts_conserved"]):
            raise ValueError(f"Incomplete or mismatched fast alignments: {name}")
        for source in SOURCES:
            run = work/source/name
            run.mkdir(parents=True)
            pair = dict(conditions=[s["arm"]], alignment_source=str(folder.parent),
                        maxi_method="manuscript", target_filter="manuscript")
            if source == "default":
                pair["published"] = str(refs/"published")
            else:
                pair["rebuilt_reference"] = rebuilt
            (run/"pair.json").write_text(json.dumps(pair))
    panel = sorted(targets("manuscript"))
    mapping = json.loads((refs/"target_guide_mapping.json").read_text())
    if len(panel) != 79 or sorted(r["target_id"] for r in mapping) != panel:
        raise ValueError("Fast and rebuilt target panels differ")
    return fast, manuscript, panel, rebuilt


def count_rows(sample, target, counts, before=None):
    if set(counts) != set(STRATEGIES):
        raise ValueError(f"Incomplete strategies for {target}")
    for strategy in STRATEGIES:
        n, d = counts[strategy]
        family_strategy = "Dedud+Dedup" if strategy.startswith("Dedud") else "Dedup"
        families = counts[family_strategy][1]
        if not (0 <= n <= d and 0 <= families <= d):
            raise ValueError(f"Invalid counts for {target}/{strategy}")
        row = dict(target=target, strategy=strategy, donor=sample["donor"], replicate=sample["replicate"],
                   arm=sample["arm"], indel=n, denominator=d, family_rows=families)
        if before is not None:
            row["pre_step6_indel"], row["pre_step6_denominator"] = before[strategy]
        yield row


def count_fast(index, target, sample, run, data):
    cmd = [sys.executable, str(AGENT/"workers/target.py"), str(index)]
    proc = subprocess.run(cmd, cwd=run, env=dict(os.environ, ABSOLVE_DATA=str(data)),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode:
        raise RuntimeError(f"Fast recount failed for {sample['sample']}/{target}: {proc.stdout[-3000:]}")
    rec = json.loads((run/sample["arm"]/"targets"/target/"counts.json").read_text())
    if (rec["target"], rec["condition"], rec["target_filter"]) != (target, sample["arm"], "manuscript"):
        raise ValueError(f"Unexpected fast count identity for {target}")
    return list(count_rows(sample, target, rec["counts"]))


def pe_families(folder, name, target):
    meta = json.loads((folder/"alignment.json").read_text())
    if (meta["sample"], meta["target"]) != (name, target):
        raise ValueError(f"Misidentified manuscript alignment: {folder}")
    fastq = folder/"CRISPResso_output.fastq.gz"
    if meta["alignment_sha256"] is None:
        if fastq.exists() or (meta["input_pairs"] and not (folder/"no-alignment.log").is_file()):
            raise ValueError(f"Inconsistent empty alignment: {folder}")
        return {}
    if sha256(fastq) != meta["alignment_sha256"]:
        raise ValueError(f"Manuscript alignment checksum differs: {folder}")
    return read_aligned_families(folder, target, meta["input_pairs"])


def recount_target(task, cfg, work, fast, manuscript, rebuilt):
    index, target = task
    refs = cfg["data"]/"references"
    identity = TargetIdentity(refs/"target_info/target_oligos_sequenes.csv")
    blocks = {}
    for source in SOURCES:
        rows = []
        for s in fast:
            rows.extend(count_fast(index, target, s, work/source/s["sample"], cfg["data"]))
        blocks["fast", source] = rows

        if source == "default":
            labels = {}
            with (refs/"published"/(target+".tsv")).open() as f:
                for barcode, flag in csv.reader(f, delimiter="\t"):
                    labels.setdefault(barcode, flag == "1")
        else:
            labels = rebuilt_labels(rebuilt, target)
        rows = []
        for s in manuscript:
            name = s["sample_name"]
            folder = cfg["manuscript_alignments"]/name/target
            families = pe_families(folder, name, target)
            counts, before = summarize(target, families, labels, identity)
            rows.extend(count_rows(s, target, counts, before))
        blocks["manuscript", source] = rows
    return blocks


def recount(cfg, work, workers):
    fast, manuscript, panel, rebuilt = prepare(cfg, work)
    run = partial(recount_target, cfg=cfg, work=work, fast=fast, manuscript=manuscript, rebuilt=rebuilt)
    merged = {(mode, source): [] for mode in ("fast", "manuscript") for source in SOURCES}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for blocks in pool.map(run, enumerate(panel)):
            for key, rows in blocks.items():
                merged[key].extend(rows)
    # The archived manuscript export visits libraries first; preserve its original byte order.
    order = {(s["donor"], s["replicate"], s["arm"]): i for i, s in enumerate(manuscript)}
    for source in SOURCES:
        merged["manuscript", source].sort(key=lambda r: order[r["donor"], r["replicate"], r["arm"]])
    return merged
