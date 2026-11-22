import gzip
import json
import os
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

import paired

from panel_report.manuscript_reads import demultiplex, summarize


class HistoricalReadTests(unittest.TestCase):
    def test_cli_routes_raw_panel_to_the_historical_membership(self):
        argv = ["paired.py", "--manuscript-panel", "--read-mode", "historical", "--fresh"]
        with patch.object(sys, "argv", argv), patch.dict(os.environ, {"ABSOLVE_HISTORICAL_MAXI_REFERENCE": "/maxi"}), \
                patch("panel_report.manuscript_run.start") as start:
            paired.main()
        args = start.call_args.args[0]
        self.assertEqual((args.maxi_method, args.maxi_reference, args.fresh), ("historical", "/maxi", True))

    def test_motif_first_occurrence_and_arm_specific_split(self):
        motif = b"CCAACCTCATAGAACA"
        const = motif + b"CTCATCC"
        tbc = b"ACGTACGTA"
        contexts = [(b"A"*16, const), (b"A"*16, motif + b"AAAAAAA"),
                    (motif, const), (b"A"*16, b"N" + const[1:])]
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "target_oligos_sequenes.csv").write_text("Barcode,OT\nACGTACGTA,focal\n")
            (root / "trim").mkdir()
            for mate in (1, 2):
                with gzip.open(root / "trim" / f"R{mate}.gz", "wb") as out:
                    for i, (pbc, fixed) in enumerate(contexts):
                        seq = pbc + fixed + tbc
                        out.write(f"@{i}/{mate}\n".encode() + seq + b"\n+\n" + b"I"*len(seq) + b"\n")
            for profile, n in (("historical-wt", 3), ("historical-hifi-noep", 2), ("manuscript", 3)):
                dest = root / profile
                dest.mkdir()
                demultiplex(dict(reference=str(root), library_profile=profile), root, dest)
                self.assertEqual(json.loads((dest / "pairs.json").read_text()), {"focal": n})
                reads = (dest / "focal.R1.fq").read_bytes().splitlines()[1::4]
                self.assertEqual(reads, [contexts[i][1] + tbc for i in (
                    (0, 1) if profile == "historical-hifi-noep" else (0, 2, 3))])

    def test_first_surviving_tie_is_used_after_cross_target_filter(self):
        discarded, edited, wt = ("R", "Z", True), ("R", "C", True), ("R", "A", False)
        families = {"pbc": {discarded: 3, edited: 2, wt: 2}}
        identity = Mock()
        identity.retain.side_effect = lambda target, a: a != discarded
        for encounter, expected in ((True, 1), (False, 0)):
            after, before = summarize("target", families, {"pbc": True}, identity, encounter)
            self.assertEqual(before["Dedud+Dedup"], [1, 1])
            self.assertEqual(after["Dedud+Dedup"], [expected, 1])
            self.assertEqual(after["Dedud"], [2, 4])
