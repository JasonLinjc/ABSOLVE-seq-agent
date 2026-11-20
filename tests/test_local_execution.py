import os
from pathlib import Path
import sys
import unittest
from unittest.mock import call, patch

import paired
from panel_report import manuscript_reads, manuscript_run


class LocalExecutionTest(unittest.TestCase):
    pair = {"NoEP": {}, "HiFi": {}}

    def invoke_pair(self, *args):
        argv = ["paired.py", "--donor", "D3734", "--replicate", "T1", *args]
        with patch.object(sys, "argv", argv), patch.object(paired, "discover",
                return_value={("D3734", "T1"): self.pair}), patch.object(paired, "execute") as execute:
            paired.main()
        return execute.call_args.args

    def test_auto_uses_local_workers_without_slurm(self):
        with patch.dict(os.environ, {}, clear=True), patch("settings.shutil.which", return_value=None):
            _, samples, workers, job = self.invoke_pair("--threads", "2")
            self.assertEqual((samples, workers, job), (self.pair, 2, None))

    def test_existing_allocation_and_local_override(self):
        with patch.dict(os.environ, {"SLURM_JOB_ID": "55", "SLURM_CPUS_PER_TASK": "8"}):
            self.assertEqual(self.invoke_pair()[2:], (8, "55"))
            self.assertEqual(self.invoke_pair("--executor", "local")[2:], (1, None))

    def test_login_node_does_not_start_fast_workers(self):
        with patch.dict(os.environ, {}, clear=True), patch("settings.shutil.which", return_value="/bin/sbatch"):
            with self.assertRaises(SystemExit):
                self.invoke_pair()

    def test_manuscript_panel_uses_same_local_stages(self):
        with patch.dict(os.environ, {}, clear=True), patch("settings.shutil.which", return_value=None), \
                patch.object(sys, "argv", ["paired.py", "--manuscript-panel", "--threads", "3"]), \
                patch.object(manuscript_run, "start") as start:
            paired.main()
        self.assertEqual(start.call_args.kwargs, {"local": True, "workers": 3})
        config = Path("/test/run.json")
        with patch.object(manuscript_reads, "worker") as worker, patch.object(manuscript_run, "finish") as finish:
            manuscript_run.run_local(config, 26, 3)
        self.assertEqual(worker.call_args_list, [
            call(config, stage, i, threads=3) for stage in ("prepare", "align") for i in range(26)])
        finish.assert_called_once_with(config)

    def test_worker_failure_stops_before_aggregation(self):
        with patch.object(manuscript_reads, "worker", side_effect=RuntimeError("failed")), \
                patch.object(manuscript_run, "finish") as finish, self.assertRaises(RuntimeError):
            manuscript_run.run_local(Path("run.json"), 2, 1)
        finish.assert_not_called()

    def test_historical_native_and_image_commands(self):
        self.assertEqual(manuscript_reads.crispresso_command({"crispresso": "/native/CRISPResso"}),
                         ["/native/CRISPResso"])
        self.assertEqual(manuscript_reads.crispresso_command({"image": "historic.sif"}, Path("/reads")),
                         ["apptainer", "exec", "--cleanenv", "--bind", "/reads", "historic.sif", "CRISPResso"])
