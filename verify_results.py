"""Check original-study outputs against the accepted Luminary and Figure 3b tables."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from zipfile import BadZipFile


# Digests come from the accepted study runs, independently of the run being tested.
PAIR = "paired_target_rates.tsv"
COUNTS = "library_counts.tsv"
PLOT = "plotted_data.tsv"
PANEL = "panel_comparison.tsv"
FIELDS = "model_field_comparison.tsv"
BASELINES = {
    "pair-fast": {PAIR: "0b0465ecc5aa244cc2a10654e153dbe8467322cbedd4f5bfc72806058018980b"},
    "pair-historical": {
        PAIR: "8cc8ca5c19e34a7838d03ee5d4658cae50b3c27664bab9959956c10faa4a22c1",
        COUNTS: "906aaa8328633745b3be80f4ff8661f4da78b046b39d4205437b6634eb1f1226",
    },
    "panel-refit": {
        COUNTS: "285a86b3b72dc352589ecf4c1df101c3bf885b3080ba14e52ed4168216714018",
        PANEL: "13a906d81fa96b325e9bdc663b8ddc3ac15c373b5861c4c9148d48de5ceaa53f",
        FIELDS: "f759f5dfd17d69a8e24b77058dafafcda91eeff5c52fcadeb3cbac3240fb1e79",
    },
    "figure3b-submitted": {PLOT: "8d3baa2324357c24f601d6075a6cc7dd45865eef166ebc47ad9a50a777652144"},
    "figure3b-refit": {PLOT: "2de36b2b73c5b59fdbbbf26c8a0bcf0999f049dc0f14a88f077d8cdf49a05c80"},
}


def check_plot_sheet(sheet, rows):
    cells = iter(sheet.values)
    headers = tuple(rows[0])
    numeric = {"order", "delta_pp", "lower_pp", "upper_pp", "q_bh", "power001", "power01"}
    if next(cells, None) != headers:
        raise ValueError("Excel source columns differ from TSV")
    for row_num, row in enumerate(rows, 2):
        values = next(cells, ())
        if len(values) != len(headers):
            raise ValueError(f"Missing or incomplete Excel row {row_num}")
        for key, cell in zip(headers, values):
            expected = row[key]
            if key in numeric and expected != "":
                matches = type(cell) in (int, float) and math.isclose(cell, float(expected), rel_tol=2e-15)
            elif key == "positive_control":
                matches = type(cell) is bool and cell == (expected == "True")
            else:
                matches = (cell if cell is not None else "") == expected
            if not matches:
                raise ValueError(f"Excel row {row_num}, {key}: expected {expected}, found {cell}")
    if next(cells, None) is not None:
        raise ValueError("Extra rows in the Figure 3b Excel table")


def check_excel(folder, rows, receipt):
    from openpyxl import load_workbook

    book = load_workbook(folder / "Figure3b_source_data.xlsx", read_only=True)
    try:
        check_plot_sheet(book["Figure 3b"], rows)
        notes = dict(list(book["Notes"].values)[1:])
        expected = {"estimates": receipt["estimates"], "power": receipt["power"],
                    "plotted_data.tsv SHA-256": receipt["plotted_data_sha256"]}
        for name, digest in receipt["input_sha256"].items():
            expected[name + " SHA-256"] = digest
        if any(notes.get(k) != v for k, v in expected.items()):
            raise ValueError("Excel source notes differ from provenance.json")
    finally:
        book.close()


def verify_power(folder):
    from analysis_report import check_coverage, read_table
    from openpyxl import load_workbook

    receipt = json.loads((folder / "power_provenance.json").read_text())
    panel = (folder / PANEL).is_file()
    sources = {PANEL, FIELDS, COUNTS} if panel else {PAIR}
    if set(receipt["input_sha256"]) != sources:
        raise ValueError("Power input inventory differs from the analysis design")
    for name, expected in receipt["input_sha256"].items():
        if hashlib.sha256((folder / name).read_bytes()).hexdigest() != expected:
            raise ValueError("Power input digest mismatch: " + name)
    rates = read_table(folder / (PANEL if panel else PAIR))
    check_coverage(rates, panel)
    effects = tuple(receipt["effects_pp"])
    if effects != (0.01, 0.1) or receipt["alpha"] != 0.05 or receipt["simulations"] <= 0:
        raise ValueError("Unexpected power simulation settings")
    keys = {(r["target"], r["strategy"], r.get("editor", "HiFi"), effect) for r in rates for effect in effects}
    rows = read_table(folder / "power_estimates.tsv")
    powers = {}
    for r in rows:
        key = r["target"], r["strategy"], r["editor"], float(r["effect_pp"])
        if r["status"] == "estimated":
            n, hits = int(r["simulations"]), int(r["detections"])
            p = float(r["power"])
            if not (n == receipt["simulations"] and 0 <= hits <= n and p == hits / n
                    and 0 <= float(r["mc_low95"]) <= p <= float(r["mc_high95"]) <= 1):
                raise ValueError(f"Invalid power estimate for {key}")
        elif (r["status"] in {"zero denominator", "effect exceeds probability range"} and r["simulations"] == "0"
              and all(r[col] == "" for col in ("power", "detections", "mc_low95", "mc_high95"))):
            p = None
        else:
            raise ValueError(f"Invalid simulation status for {key}")
        powers[key] = p
    if set(powers) != keys or len(rows) != len(keys):
        raise ValueError("Power table lacks complete unique coverage")

    figures = folder / "figures"
    plotted = read_table(figures / PLOT)
    deltas = {(r["target"], r.get("editor", "HiFi")): r["delta_pp"] if panel else r["delta_HiFi_minus_NoEP_pp"]
              for r in rates if r["strategy"] == "Dedud"}
    if len(plotted) != len(deltas) or {(r["target"], r["editor"]) for r in plotted} != set(deltas):
        raise ValueError("Incomplete plotted target/editor records")
    for r in plotted:
        if r["strategy"] != "Dedud":
            raise ValueError("Unexpected plotted strategy")
        expected = [deltas[r["target"], r["editor"]]] + [powers[r["target"], "Dedud", r["editor"], e] for e in effects]
        for col, val in zip(("delta_pp", "power001", "power01"), expected):
            actual = r[col]
            if val is None or val == "":
                matches = actual == ""
            else:
                matches = actual != "" and math.isclose(float(actual), float(val), rel_tol=1e-12, abs_tol=1e-12)
            if not matches:
                raise ValueError("Plotted value differs from run tables: " + col)
    stem = "Figure3b" if panel else "paired_effect_power"
    workbook = figures / (("Figure3b_source_data" if panel else stem) + ".xlsx")
    book = load_workbook(workbook, read_only=True, data_only=True)
    try:
        check_plot_sheet(book.worksheets[0], plotted)
    finally:
        book.close()
    if not (figures / (stem + ".pdf")).read_bytes().startswith(b"%PDF-") or not (
            figures / (stem + ".png")).read_bytes().startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("Missing or invalid PDF/PNG figure")


def verify(profile, folder):
    if profile == "power":
        return verify_power(folder)
    for filename, expected in BASELINES[profile].items():
        actual = hashlib.sha256((folder / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"{filename}: expected SHA-256 {expected}, found {actual}")
    if profile.startswith("figure3b-"):
        refit = profile == "figure3b-refit"
        sources = {
            "SupplementaryTable2.xlsx": "c92271bab3d95eb729fc6900de1d2024d12c13a4622abc01d0aa934d6b19f4af",
            "target_order.tsv": "dc5a1b056d979e7dc4913d28e05286c8843b86b8cb786f8c89de2c2593bdbe57",
        }
        if refit:
            sources.update({k: BASELINES["panel-refit"][k] for k in (PANEL, FIELDS)})
        receipt = json.loads((folder / "provenance.json").read_text())
        expected = dict(figure="3b", targets=79, rows=158, input_sha256=sources,
                        estimates="agent refit" if refit else "submitted workbook",
                        power="submitted candidate-specific simulations", plotted_data_sha256=BASELINES[profile][PLOT])
        if any(receipt.get(k) != v for k, v in expected.items()):
            raise ValueError("Figure 3b provenance differs from the selected benchmark")
        with (folder / PLOT).open(newline="") as stream:
            check_excel(folder, list(csv.DictReader(stream, delimiter="\t")), receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=[*BASELINES, "power", "panel-historical"])
    parser.add_argument("directory", type=Path, help="Directory containing the exported tables")
    parser.add_argument("--workbooks", type=Path, help="Submitted workbook directory for panel-historical")
    args = parser.parse_args()
    if (args.profile == "panel-historical") != bool(args.workbooks):
        parser.error("--workbooks is required for panel-historical and only applies to that profile")
    try:
        if args.profile == "panel-historical":
            from panel_report.check_historical_panel import verify as verify_historical
            verify_historical(args.directory, args.workbooks)
        else:
            verify(args.profile, args.directory)
    except (OSError, ValueError, KeyError, BadZipFile) as error:
        parser.exit(1, f"FAIL {args.profile}: {error}\n")
    print(f"PASS {args.profile}")


if __name__ == "__main__":
    main()
