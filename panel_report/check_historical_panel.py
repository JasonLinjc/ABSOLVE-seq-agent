"""Check a fresh historical panel against the four submitted workbooks."""
from collections import Counter
import hashlib
import json
import math

from openpyxl import load_workbook

from analysis_report import check_coverage, read_table


def bh(values):
    ordered = sorted(values, key=values.get)
    adjusted = {}
    upper = 1.0
    for rank in range(len(ordered), 0, -1):
        target = ordered[rank - 1]
        upper = min(upper, values[target] * len(ordered) / rank)
        adjusted[target] = upper
    return adjusted


def verify(folder, workbooks):
    receipt = json.loads((folder / "provenance.json").read_text())
    if (receipt["read_mode"], receipt["maxi_method"], receipt["libraries"], receipt["targets"]) != (
            "historical", "historical", 26, 79):
        raise ValueError("Expected a historical raw-read panel")
    inventory = read_table(folder / "sample_inventory.tsv")
    libraries = {(r["arm"], r["donor"], r["replicate"]) for r in inventory}
    if len(libraries) != 26 or Counter(r["arm"] for r in inventory) != {"WT": 9, "HiFi": 8, "NoEP": 9}:
        raise ValueError("Incomplete historical library inventory")
    rows = read_table(folder / "panel_comparison.tsv")
    check_coverage(rows, panel=True)
    results = {(r["strategy"], r["editor"], r["target"]): r for r in rows}
    counts = read_table(folder / "library_counts.tsv")
    expected = {(s, t, *lib) for s, e, t in results for lib in libraries}
    actual = {(r["strategy"], r["target"], r["arm"], r["donor"], r["replicate"]) for r in counts}
    if actual != expected or len(actual) != len(counts):
        raise ValueError("Incomplete or duplicate library counts")
    matches = positives = 0
    errors = dict(control_pp=0.0, edited_pp=0.0, delta_pp=0.0, q=0.0)
    for num, strategy in enumerate(("Dedud", "Unfiltered", "Dedup", "Dedud+Dedup"), 2):
        path = workbooks / f"SupplementaryTable{num}.xlsx"
        if hashlib.sha256(path.read_bytes()).hexdigest() != receipt["submitted_workbooks"][path.name]:
            raise ValueError(f"Workbook differs from the model reference: {path.name}")
        wb = load_workbook(path, read_only=True, data_only=True)
        try:
            data = iter(wb.active.values)
            next(data)
            headers = next(data)
            submitted = [dict(zip(headers, row)) for row in data if row and row[0]]
        finally:
            wb.close()
        targets = {r["barcode"] for r in submitted}
        if len(submitted) != 79 or targets != {t for s, e, t in results if s == strategy}:
            raise ValueError(f"Different target panel in {path.name}")
        for editor, arm in (("WT", "WTCas9"), ("HiFi", "HiFi")):
            reference_q = bh({r["barcode"]: r["asymPval-treat" + arm] for r in submitted})
            replay_q = bh({t: float(results[strategy, editor, t]["p_one_sided"]) for t in targets})
            for ref in submitted:
                t = ref["barcode"]
                row = results[strategy, editor, t]
                control, edited = (100 * ref["predicted-" + a] for a in ("NoEP", arm))
                ref_delta = edited - control
                delta, q = float(row["delta_pp"]), replay_q[t]
                check = {"q_bh": q, "manuscript_q_bh": reference_q[t], "manuscript_delta_pp": ref_delta,
                         "delta_pp": float(row["edited_rate_pct"]) - float(row["control_rate_pct"])}
                if not all(math.isclose(float(row[k]), v, rel_tol=1e-10, abs_tol=5e-12) for k, v in check.items()):
                    raise ValueError(f"Inconsistent exported delta or q: {strategy}/{editor}/{t}")
                hit = delta > 0.1 and q < 0.05
                ref_hit = ref_delta > 0.1 and reference_q[t] < 0.05
                if (row["agent_hit"], row["submitted_hit"]) != tuple(str(val).upper() for val in (hit, ref_hit)):
                    raise ValueError(f"Inconsistent exported call: {strategy}/{editor}/{t}")
                matches += hit == ref_hit
                positives += hit
                residuals = dict(control_pp=float(row["control_rate_pct"]) - control,
                                 edited_pp=float(row["edited_rate_pct"]) - edited, delta_pp=delta - ref_delta,
                                 q=q - reference_q[t])
                for field, residual in residuals.items():
                    if not math.isfinite(residual):
                        raise ValueError(f"Nonfinite {field} for {strategy}/{editor}/{t}")
                    errors[field] = max(errors[field], abs(residual))
    print(f"{len(counts)} library counts, {matches}/632 calls agree, {positives} positive")
    print("Maximum absolute differences", json.dumps(errors))
    if matches != 632 or max(errors.values()) > 0.01:
        raise ValueError("Historical panel exceeds the accepted call or numerical tolerance")
