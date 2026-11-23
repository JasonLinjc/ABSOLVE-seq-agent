import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from examples.d3734_t1_manifest import make_manifest
from panel_report.manuscript_run import panel_references


class PanelInputTest(unittest.TestCase):
    def test_public_manifest_requires_both_archive_checksums(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            files = [root / f"SRR1_{mate}.fastq.gz" for mate in (1, 2)]
            for file in files:
                file.write_bytes(file.name.encode())
            with (root / "study_panel_ena.tsv").open("w") as catalog:
                catalog.write("sample_name\tlibrary_id\trun_accession\tr1_url\tmd5_r1\tr2_url\tmd5_r2\n")
                columns = ["WT_pool_D3742_T1", "LV_S1", "SRR1"]
                for file in files:
                    columns.extend([f"https://example.org/{file.name}", hashlib.md5(file.read_bytes()).hexdigest()])
                catalog.write("\t".join(columns) + "\n")
            with patch("examples.d3734_t1_manifest.__file__", str(root / "make_manifest.py")):
                output = root / "manifest.tsv"
                make_manifest(root, output, ena=True, panel=True)
                with output.open() as stream:
                    row, = csv.DictReader(stream, delimiter="\t")
                for mate, file in enumerate(files, 1):
                    self.assertEqual(row[f"sha256_r{mate}"], hashlib.sha256(file.read_bytes()).hexdigest())
                files[1].write_bytes(b"incomplete transfer")
                with self.assertRaisesRegex(ValueError, "ENA MD5 mismatch"):
                    make_manifest(root, root / "invalid.tsv", ena=True, panel=True)
                self.assertFalse((root / "invalid.tsv").exists())

    def test_missing_optional_and_incomplete_references(self):
        with tempfile.TemporaryDirectory() as tmp, patch("panel_report.manuscript_run.DATA", Path(tmp)):
            self.assertIsNone(panel_references())
            with self.assertRaises(FileNotFoundError):
                panel_references(required=True)
            refs = Path(tmp) / "submitted_model_reference"
            refs.mkdir()
            with self.assertRaises(FileNotFoundError):
                panel_references()
            for n in range(2, 6):
                (refs / f"SupplementaryTable{n}.csv").touch()
            (refs / "source_sha256.json").touch()
            self.assertEqual(panel_references(required=True), refs.resolve())

    def test_panel_manifest_from_original_laboratory_filenames(self):
        catalog = Path(__file__).resolve().parents[1] / "examples/study_panel_ena.tsv"
        with tempfile.TemporaryDirectory() as tmp, catalog.open() as handle:
            root = Path(tmp)
            records = list(csv.DictReader(handle, delimiter="\t"))
            for row in records:
                condition, _, donor, rep = row["sample_name"].split("_")
                folder = root / condition
                folder.mkdir(exist_ok=True)
                for mate in (1, 2):
                    if condition == "WT" and rep != "T3":
                        name = f"Nola_EP031022_2_LVpool_wtCas9_{donor}_{rep}_R{mate}.fastq.gz"
                    else:
                        name = f"{row['library_id']}_CKDL123_{mate}.fq.gz"
                    (folder / name).touch()
            out = root / "manifest.tsv"
            make_manifest(root, out, panel=True)
            with out.open() as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(rows), 26)
            self.assertEqual({r["sample_name"] for r in rows}, {r["sample_name"] for r in records})
            for r in rows:
                for mate in (1, 2):
                    self.assertEqual(root / r["pair_key"].format(mate=mate), Path(r[f"r{mate}"]))
