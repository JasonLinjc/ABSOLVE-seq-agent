import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from tools import luminary_pbc as launcher


class LuminaryDeliveryTests(unittest.TestCase):
    def test_resume_delivers_completed_differences_without_submitting(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            receipt = dict(host="luminary", job="42", remote_run="/archive/job 42", log="/archive/log")
            (dest / "job.json").write_text(json.dumps(receipt))
            files = {"summary.tsv": b"one comparison\n", "differences.tsv": b"one difference\n"}
            manifest = "".join(
                f"{hashlib.sha256(data).hexdigest()}  {name}\n" for name, data in files.items()
            ).encode()
            responses = [b"42|FAILED|1:0\n", b'{"status":"different","differences":1}',
                         manifest, *files.values()]
            with patch.object(sys, "argv", ["launcher", "--resume", tmp]), \
                 patch.object(launcher, "remote", side_effect=responses) as ssh, \
                 patch.object(launcher, "submit") as submit:
                self.assertEqual(launcher.main(), 1)
                submit.assert_not_called()
                self.assertIn("sha256sum --check", ssh.call_args_list[2].args[1])
            for name, data in files.items():
                self.assertEqual((dest / name).read_bytes(), data)
            self.assertEqual((dest / "SHA256SUMS").read_bytes(), manifest)

    def test_bad_download_does_not_publish_partial_reports(self):
        receipt = dict(host="luminary", remote_run="/archive", log="/archive/log")
        manifest = f"{hashlib.sha256(b'ok').hexdigest()}  summary.tsv\n{'0' * 64}  differences.tsv\n"
        responses = [b'{"status":"pass","differences":0}', manifest.encode(), b"ok", b"corrupt"]
        with tempfile.TemporaryDirectory() as tmp, patch.object(launcher, "remote", side_effect=responses):
            with self.assertRaisesRegex(ValueError, "Downloaded digest mismatch for differences.tsv"):
                launcher.collect(receipt, Path(tmp), "COMPLETED", "0:0")
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_job_failure_is_not_a_scientific_difference(self):
        receipt = dict(host="luminary", remote_run="/archive", log="/archive/slurm-42.out")
        with patch.object(launcher, "remote") as ssh:
            with self.assertRaisesRegex(ValueError, "slurm-42.out"):
                launcher.collect(receipt, Path("."), "FAILED", "2:0")
            ssh.assert_not_called()
