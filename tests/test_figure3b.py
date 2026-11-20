import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from openpyxl import Workbook, load_workbook

from verify_results import check_excel, verify


SKILL = Path(__file__).resolve().parents[1] / ".agents/skills/absolve-seq-figures"
SCRIPT = SKILL / "scripts/manuscript_figures.py"


class Figure3bTest(unittest.TestCase):
    def test_render_and_verify_excel(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            out = folder / "render"
            book = Workbook()
            sheet = book.active
            sheet.append(["Synthetic Figure 3b test"])
            headings = ["barcode", "predicted-NoEP"]
            for suffix in ("WTCas9", "HiFi"):
                headings.extend([name + suffix for name in ("predicted-", "confLow-", "confHigh-",
                                "asymPval-treat", "power0.01Treat", "power0.1Treat")])
            sheet.append(headings)
            with (SKILL / "references/target_order.tsv").open() as stream:
                order = list(csv.DictReader(stream, delimiter="\t"))
            for target in reversed(order):
                sheet.append([target["target"], .01, .02, .009, .04, .001, .6, .9,
                              .014, .012, .02, .15, .4, .8])
            book.save(folder / "SupplementaryTable2.xlsx")

            subprocess.run([sys.executable, str(SCRIPT), "3b", "--tables", str(folder),
                            "--out", str(out)], check=True)
            receipt = json.loads((out / "provenance.json").read_text())
            self.assertEqual(receipt["plotted_data_sha256"],
                             hashlib.sha256((out / "plotted_data.tsv").read_bytes()).hexdigest())
            self.assertEqual({p.name for p in out.iterdir()}, {"Figure3b.pdf", "Figure3b.png",
                             "Figure3b_source_data.xlsx", "provenance.json", "plotted_data.tsv"})
            with (out / "plotted_data.tsv").open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([(r["target"], r["editor"]) for r in rows],
                             [(t["target"], e) for e in ("WT", "HiFi") for t in order])
            for i, delta, low, q in [(0, 1., -.1, .001), (79, .4, .2, .15)]:
                self.assertAlmostEqual(float(rows[i]["delta_pp"]), delta)
                self.assertAlmostEqual(float(rows[i]["lower_pp"]), low)
                self.assertAlmostEqual(float(rows[i]["q_bh"]), q)
            check_excel(out, rows, receipt)
            # A valid table from a different input must fail the original-study benchmark.
            with self.assertRaisesRegex(ValueError, "expected SHA-256"):
                verify("figure3b-submitted", out)

            path = out / "Figure3b_source_data.xlsx"
            altered = load_workbook(path)
            altered["Figure 3b"]["G2"] = 0.5
            altered.save(path)
            with self.assertRaisesRegex(ValueError, "Excel row 2, delta_pp"):
                check_excel(out, rows, receipt)
