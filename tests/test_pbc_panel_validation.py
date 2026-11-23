"""Exercise the panel acceptance diff on scientific changes and serialization changes."""
import tempfile
import unittest
from pathlib import Path

from panel_report.manuscript_reads import save
from tools.validate_pbc_panel import compare


class PanelValidationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.left = Path(self.tmp.name)/"default.tsv"
        self.right = Path(self.tmp.name)/"rebuilt.tsv"
        self.rows = [dict(target=str(t), strategy=s, editor=e, q_bh="0.01")
                     for t in range(79) for s in ("Unfiltered", "Dedup", "Dedud", "Dedud+Dedup")
                     for e in ("WT", "HiFi")]
        save(self.left, self.rows)

    def check_right(self, rows):
        save(self.right, rows)
        differences = []
        summary = compare("manuscript", "default_vs_rebuilt", "equal-donor",
                          self.left, self.right, differences)
        self.assertEqual(summary["differences"], len(differences))
        return differences

    def test_identical_exports(self):
        self.assertEqual(self.check_right(self.rows), [])

    def test_changed_q_names_the_contrast(self):
        self.rows[0]["q_bh"] = "0.02"
        diff, = self.check_right(self.rows)
        self.assertEqual((diff["target"], diff["strategy"], diff["editor"], diff["field"],
                          diff["left"], diff["right"]), ("0", "Unfiltered", "WT", "q_bh", "0.01", "0.02"))

    def test_different_bytes_even_when_keys_match(self):
        diff, = self.check_right(list(reversed(self.rows)))
        self.assertEqual(diff["field"], "byte_sha256")

    def test_missing_and_duplicate_keys(self):
        self.assertEqual({d["field"] for d in self.check_right(self.rows[:-1])},
                         {"rows_expected_632", "row_present"})
        self.right.unlink()
        with self.assertRaisesRegex(ValueError, "Duplicate key"):
            self.check_right(self.rows + self.rows[:1])
