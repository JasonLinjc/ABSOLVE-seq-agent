"""Audit the two historical allele tables and both filtered counting strategies by pBC."""
import hashlib
import json

import pandas as pd

from .collect_historical_families import WINDOW, FULL

RUN_DIRS = ("singleUMI_crispresso", "singleUMI_crispresso_positive",
            "singleUMI_crispresso_positive_borderline")
METRICS = ("dedud_indel_reads", "dedud_total_reads", "dedud_dedup_indel_families",
           "dedud_dedup_total_families")
PRESENCE = ("expected", "in_replay", "in_window", "in_full")
ID_PATTERN = "[ACGTN]{16}"


def digest_ids(ids):
    return hashlib.sha256("".join(b + "\n" for b in sorted(ids)).encode()).hexdigest()


def unique_ids(values, path):
    ids = pd.Index(values, dtype=object)
    if not ids.is_unique or not ids.to_series().str.fullmatch(ID_PATTERN).fillna(False).all():
        raise ValueError(f"Invalid or duplicate pBC in {path}")
    return set(ids)


def historical_ids(root, names):
    """Accept a list only when its successes and failures account for the batch plan."""
    expected = None
    records = []
    for name in names:
        for directory in RUN_DIRS:
            batch = root / directory / ("CRISPResso_on_" + name) / ("CRISPRessoBatch_on_" + name + "_input")
            info = batch / "CRISPResso2Batch_info.json"
            if not info.is_file():
                continue
            contents = info.read_bytes()
            results = json.loads(contents).get("results", {})
            record = dict(path=str(info), sha256=hashlib.sha256(contents).hexdigest(), status="incomplete")
            records.append(record)
            plan = root / "singleUMI_crispresso_input" / (name + "_input.tsv")
            if not plan.is_file() or not {"completed_batch_arr", "failed_batch_arr"} <= results.keys():
                continue
            planned = unique_ids(pd.read_csv(plan, sep="\t", usecols=["name"]).name, plan)
            done = unique_ids(results["completed_batch_arr"], info)
            failed = unique_ids(results["failed_batch_arr"], info)
            if done & failed:
                raise ValueError(f"Overlapping successful and failed families in {info}")
            if done | failed != planned:
                record["status"] = "different_input_universe"
                continue
            quant = batch / "CRISPRessoBatch_quantification_of_editing_frequency.txt"
            q = pd.read_csv(quant, sep="\t", usecols=["Batch", "Amplicon"])
            if unique_ids(q.Batch, quant) != done or not q.Amplicon.eq("Reference").all():
                raise ValueError(f"Quantification disagrees with completed families in {info}")
            if expected is not None and done != expected:
                raise ValueError(f"Conflicting complete historical lists for {name}")
            expected = done
            record.update(status="complete", plan=str(plan), quantification=str(quant),
                          successful=len(done), successful_ids_sha256=digest_ids(done))
    return expected, records


def historical_tables(root, arm, donor, replicate, target, lv):
    token = target.replace("_", "-")
    group = "WTCas9" if arm == "WT" else arm
    raw = "singleUMI_crispresso_WTCas9" if arm == "WT" else "alleles_all_HiFi_NoEP"
    prefix = f"{token}_{lv.replace('_', '-')}_{donor}_{group}_{replicate}"
    paths = [root / raw / "alleles_all" / (prefix + suffix) for suffix in (WINDOW, FULL)]
    window, full = (set(pd.read_csv(p, sep="\t", usecols=["UMI"]).UMI) for p in paths)
    filename = f"{target}_{lv}_{donor}_{group}_{replicate}_Allele_frequency_table_full_withUMI.annot.txt"
    counts = []
    for directory in ("dedudS_filtered_v3", "dedudS_filtered_v3_dedup"):
        path = root / "target_umi_barcode_table_WTCas9_HiFi_NoEP" / directory / filename
        paths.append(path)
        frame = pd.read_csv(path, sep="\t", usecols=["UMI", "#Reads", "Unedited"],
                            dtype={"#Reads": "int64", "Unedited": "boolean"})
        if frame.isna().any().any() or (frame["#Reads"] < 1).any():
            raise ValueError(f"Invalid filtered allele in {path}")
        if directory.endswith("_dedup"):
            if not frame.UMI.is_unique:
                raise ValueError(f"Multiple representatives per pBC in {path}")
            frame["#Reads"] = 1
        frame["indel"] = frame["#Reads"] * (~frame.Unedited).astype("int64")
        counts.append(frame.groupby("UMI")[["indel", "#Reads"]].sum())
    archive = pd.concat(counts, axis=1)
    archive.columns = list(METRICS)
    return window, full, archive, paths


def replay_counts(base, arm, donor, replicate, target, reference_hash):
    table = base / (target + ".tsv.gz")
    receipt_path = base / (target + ".json")
    receipt = json.loads(receipt_path.read_bytes())
    expected_meta = dict(target=target, donor=donor, replicate=replicate,
                         read_profile="historical-wt" if arm == "WT" else "historical-hifi-noep",
                         maxi_reference_sha256=reference_hash)
    if receipt.get("arm", arm) != arm or any(receipt.get(k) != v for k, v in expected_meta.items()):
        raise ValueError(f"Replay metadata do not match the requested batch in {receipt_path}")
    f = pd.read_csv(table, sep="\t", dtype={"aligned": "int64", "indel": "int64", "eligible": "int64",
                    "dedud_reads": "int64", "dedud_indel": "int64", "retained_representative_indel": "Int64"}
                    ).set_index("pbc", verify_integrity=True)
    ids = unique_ids(f.index, table)
    retained = f.dedud_reads.gt(0)
    representative = f.retained_representative_indel
    if (not f.eligible.isin([0, 1]).all() or (f.aligned < 1).any() or
            not f.indel.between(0, f.aligned).all() or not f.dedud_reads.between(0, f.aligned).all() or
            not f.dedud_indel.between(0, f.dedud_reads).all() or (f.dedud_indel > f.indel).any() or
            (retained & f.eligible.eq(0)).any() or not representative.notna().equals(retained) or
            not representative.dropna().isin([0, 1]).all()):
        raise ValueError(f"Invalid family counts in {table}")
    counts = receipt["counts"]
    if (receipt["families"] != len(ids) or digest_ids(ids) != receipt["sha256"]
            or len(ids) != counts["Dedup"][1] or [f.indel.sum(), f.aligned.sum()] != counts["Unfiltered"]
            or [f.dedud_indel.sum(), f.dedud_reads.sum()] != counts["Dedud"]
            or [representative.sum(), retained.sum()] != counts["Dedud+Dedup"]):
        raise ValueError(f"Replay family IDs or counts disagree with {receipt_path}")
    frame = pd.concat([f.dedud_indel, f.dedud_reads, representative.fillna(0), retained], axis=1)
    frame.columns = list(METRICS)
    return ids, frame.astype("int64"), table, receipt_path


def compare(expected, rebuilt, window, full, archive, replay):
    index = sorted(expected | rebuilt | window | full | set(archive.index))
    unique_ids(index, "family comparison")
    rows = pd.DataFrame(index=pd.Index(index, name="pbc"))
    for name, ids in zip(PRESENCE, [expected, rebuilt, window, full]):
        rows[name] = rows.index.isin(ids)
    flag = ((rows.expected != rows.in_replay) | (rows.expected != rows.in_window)
            | (rows.expected != rows.in_full))
    summary = dict(expected_families=len(expected), replay_families=len(rebuilt),
                   window_families=len(window), full_families=len(full),
                   missing_window=len(expected - window), missing_full=len(expected - full),
                   unexpected_window=len(window - expected), unexpected_full=len(full - expected),
                   missing_replay=len(expected - rebuilt), unexpected_replay=len(rebuilt - expected))
    archived = archive.reindex(index).fillna(0).astype("int64")
    current = replay.reindex(index).fillna(0).astype("int64")
    for metric in METRICS:
        delta = current[metric] - archived[metric]
        for label, values in (("archived", archived[metric]), ("replay", current[metric]), ("delta", delta)):
            column = label + "_" + metric
            rows[column] = values
            summary[column] = int(values.sum())
        flag |= delta.ne(0)
    return summary, rows.loc[flag].reset_index()


def audit_batch(task):
    root, base, arm, donor, replicate, target, lv, ref_hash = task
    key = dict(arm=arm, donor=donor, replicate=replicate, target=target)
    try:
        token = target.replace("_", "-")
        names = [f"{arm}_{donor}_{replicate}_{token}", f"{lv.replace('_', '-')}_{token}"]
        direct, records = historical_ids(root, names)
        rebuilt, replay, table, receipt = replay_counts(base, arm, donor, replicate, target, ref_hash)
        window, full, archive, paths = historical_tables(root, arm, donor, replicate, target, lv)
        expected = rebuilt if direct is None else direct
        key["evidence"] = "fastq_reconstruction" if direct is None else "historical_success_list"
        summary, differences = compare(expected, rebuilt, window, full, archive, replay)
        source = dict(**key, expected_ids_sha256=digest_ids(expected), historical_records=records,
                      replay_table=str(table), replay_receipt=str(receipt),
                      archived_tables=dict(zip(("window", "full", "dedud", "dedud_dedup"), map(str, paths))))
        return dict(**key, **summary, discrepant_pbcs=len(differences)), [
            dict(**key, **r) for r in differences.to_dict("records")], source
    except Exception as error:
        raise ValueError(f"{arm}/{donor}/{replicate}/{target}: {error}") from error
