import csv
import json
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest

from analysis_report import STRATEGIES, write_report


def write_table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class RunReportTest(unittest.TestCase):
    targets = ["1617_OT_0040_ALT"] + [f"target_{i}" for i in range(78)]

    def test_paired_report_counts_zero_denominators_and_checks_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            rows = []
            for target in self.targets:
                observed = target == self.targets[0]
                for strategy in STRATEGIES:
                    rows.append(dict(target=target, strategy=strategy, donor="D1", replicate="T1",
                                     maxi_method="manuscript", target_filter="manuscript",
                                     NoEP_indel=0, NoEP_denominator=10 if observed else 0,
                                     NoEP_rate_pct="0" if observed else "",
                                     HiFi_indel=int(observed), HiFi_denominator=10 if observed else 0,
                                     HiFi_rate_pct="10" if observed else "",
                                     delta_HiFi_minus_NoEP_pp="10" if observed else ""))
            table = dest / "paired_target_rates.tsv"
            write_table(table, rows)
            qc = {arm: dict(sample_id=arm, pairs=50, prepared_pairs=25) for arm in ("HiFi", "NoEP")}
            (dest / "qc.json").write_text(json.dumps(qc))
            write_report(dest, "fast")
            report = (dest / "report.md").read_text()
            self.assertIn("| NoEP | NoEP | 50 | 25 | 50 |", report)
            self.assertIn("| Dedud | 0 | 10 | 10 | 78 |", report)
            self.assertIn("| Dedud | 0 / 10 | 0 | 1 / 10 | 10 | 10 |", report)
            self.assertIn("[qc.json](qc.json)", report)
            write_table(table, rows[:-1] + [rows[0]])
            with self.assertRaisesRegex(ValueError, "duplicated"):
                write_report(dest, "fast")

    def test_panel_report_displays_discordances_and_input_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            rows = [dict(target=t, strategy=s, editor=e, agent_hit="FALSE", submitted_hit="FALSE",
                         delta_pp="0", manuscript_delta_pp="0", q_bh="1", manuscript_q_bh="1")
                    for t in self.targets for s in STRATEGIES for e in ("WT", "HiFi")]
            rows[0].update(agent_hit="TRUE", delta_pp="0.2", q_bh="0.01")
            write_table(dest / "panel_comparison.tsv", rows)
            write_report(dest, "frozen-count refit", panel=True)
            report = (dest / "report.md").read_text()
            self.assertIn("frozen-count refit", report)
            self.assertIn("631 of 632", report)
            self.assertIn("0.2 pp for Δ and 0.99 for q", report)
            self.assertIn("| Unfiltered | WT | 1 | 0 | 1 |", report)
            self.assertIn("| 1617_OT_0040_ALT | Unfiltered | WT | 0.2 | 0 | 0.01 | 1 |", report)

    def test_fast_workflow_publishes_tables_and_report_together(self):
        script = Path(__file__).resolve().parents[1] / "workers/summarize.py"
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            config = dict(donor="D1", replicate="T1", maxi_method="manuscript", min_exact_reads=None,
                          target_filter="manuscript")
            (run / "pair.json").write_text(json.dumps(config))
            for arm in ("NoEP", "HiFi"):
                prepared = run / arm / "prepared"
                prepared.mkdir(parents=True)
                qc = dict(status="complete", scope="full", counts_conserved=True, sample_id=arm,
                          pairs=79, prepared_pairs=79, dispositions={}, structures={},
                          input_config_sha256={}, source_sha256={})
                (prepared / "summary.json").write_text(json.dumps(qc))
                for target in self.targets:
                    folder = run / arm / "targets" / target
                    folder.mkdir(parents=True)
                    counts = dict(target=target, condition=arm, prepared=1,
                                  counts={s: [0, 1] for s in STRATEGIES}, mixed_top_ties={})
                    (folder / "counts.json").write_text(json.dumps(counts))
            process = subprocess.run([sys.executable, str(script), "0"], cwd=run, capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            report = (run / "results/report.md").read_text()
            self.assertIn("D1 / T1", report)
            self.assertIn("[Run configuration](../pair.json)", report)
            self.assertTrue((run / "results/paired_target_rates.tsv").is_file())
