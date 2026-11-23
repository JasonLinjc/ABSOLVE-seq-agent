import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from panel_report.family_audit import METRICS, RUN_DIRS, compare, historical_ids, replay_counts


class FamilyAuditTest(unittest.TestCase):
    def test_per_table_membership_and_offsetting_counts(self):
        a, b, c = ("A" * 16, "C" * 16, "G" * 16)
        archived = pd.DataFrame([[1, 5, 1, 1], [0, 4, 0, 1]], index=[a, b], columns=METRICS)
        replay = pd.DataFrame([[0, 4, 0, 1], [1, 5, 1, 1]], index=[a, b], columns=METRICS)
        row, detail = compare({a, b}, {a, b}, {a}, {a, b, c}, archived, replay)
        self.assertEqual((row["missing_window"], row["unexpected_full"]), (1, 1))
        self.assertEqual(row["delta_dedud_total_reads"], 0)
        self.assertEqual(detail.set_index("pbc").delta_dedud_total_reads.to_dict(), {a: -1, b: 1, c: 0})
        self.assertFalse(detail.set_index("pbc").loc[b, "in_window"])

    def test_incomplete_list_and_conflicting_complete_sources(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, b = "A" * 16, "C" * 16
            plan = root / "singleUMI_crispresso_input" / "sample_input.tsv"
            plan.parent.mkdir()
            plan.write_text(f"name\n{a}\n{b}\n")

            def source(directory, done, failed):
                batch = root / directory / "CRISPResso_on_sample" / "CRISPRessoBatch_on_sample_input"
                batch.mkdir(parents=True, exist_ok=True)
                results = dict(completed_batch_arr=done, failed_batch_arr=failed)
                (batch / "CRISPResso2Batch_info.json").write_text(json.dumps(dict(results=results)))
                (batch / "CRISPRessoBatch_quantification_of_editing_frequency.txt").write_text(
                    "Batch\tAmplicon\n" + "".join(f"{v}\tReference\n" for v in done))

            source(RUN_DIRS[0], [a], [])
            expected, records = historical_ids(root, ["sample"])
            self.assertIsNone(expected)
            self.assertEqual(records[0]["status"], "different_input_universe")
            source(RUN_DIRS[1], [a], [b])
            self.assertEqual(historical_ids(root, ["sample"])[0], {a})
            source(RUN_DIRS[2], [b], [a])
            with self.assertRaisesRegex(ValueError, "Conflicting"):
                historical_ids(root, ["sample"])

    def test_receipt_rejects_swapped_library(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "target.json").write_text(json.dumps({"donor": "D3908"}))
            with self.assertRaisesRegex(ValueError, "metadata"):
                replay_counts(root, "WT", "D3734", "T1", "target", "hash")
