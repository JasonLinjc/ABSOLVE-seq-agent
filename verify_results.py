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


def check_excel(folder, rows, receipt):
    from openpyxl import load_workbook

    book = load_workbook(folder / "Figure3b_source_data.xlsx", read_only=True)
    numeric = {"delta_pp", "lower_pp", "upper_pp", "q_bh", "power001", "power01"}
    try:
        data = iter(book["Figure 3b"].values)
        headers = tuple(rows[0])
        if next(data, None) != headers:
            raise ValueError("Figure 3b Excel columns differ from the TSV")
        for row_num, row in enumerate(rows, 2):
            cells = next(data, ())
            if len(cells) != len(headers):
                raise ValueError(f"Missing or incomplete Excel row {row_num}")
            for key, cell in zip(headers, cells):
                value = row[key]
                if key in numeric:
                    matches = type(cell) in (float, int) and math.isclose(cell, float(value), rel_tol=2e-15)
                elif key == "positive_control":
                    matches = type(cell) is bool and cell == (value == "True")
                elif key == "order":
                    matches = type(cell) is int and cell == int(value)
                else:
                    matches = cell == value
                if not matches:
                    raise ValueError(f"Excel row {row_num}, {key}: expected {value}, found {cell}")
        if next(data, None) is not None:
            raise ValueError("Extra rows in the Figure 3b Excel table")
        notes = dict(list(book["Notes"].values)[1:])
        expected = {"estimates": receipt["estimates"], "power": receipt["power"],
                    "plotted_data.tsv SHA-256": receipt["plotted_data_sha256"]}
        for name, digest in receipt["input_sha256"].items():
            expected[name + " SHA-256"] = digest
        if any(notes.get(k) != v for k, v in expected.items()):
            raise ValueError("Excel source notes differ from provenance.json")
    finally:
        book.close()


def verify(profile, folder):
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
    parser.add_argument("profile", choices=BASELINES)
    parser.add_argument("directory", type=Path, help="Directory containing the exported tables")
    args = parser.parse_args()
    try:
        verify(args.profile, args.directory)
    except (OSError, ValueError, KeyError, BadZipFile) as error:
        parser.exit(1, f"FAIL {args.profile}: {error}\n")
    print(f"PASS {args.profile}")


if __name__ == "__main__":
    main()
