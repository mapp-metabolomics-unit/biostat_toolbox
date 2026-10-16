"""Metadata matching and end-to-end mzmlwatcher enrichment checks."""

import csv
import hashlib
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/enrich_metadata.py"
FIXTURE = Path(__file__).parent / "fixtures/enrichment.mzML"
spec = importlib.util.spec_from_file_location("enrich_metadata", SCRIPT)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
META_COLUMNS = ["filename", "sample_id", "sample_type"]
SCAN_COLUMNS = [
    "file_path", "file_name", "parse_status", "acquisition_date", "polarity",
    "raw_file_sha1", "converted_file_sha1", "checksum",
]


def write_table(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class EnrichMetadataTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        root = Path(self.folder.name)
        self.metadata = root / "metadata.tsv"
        self.watcher = root / "watcher.tsv"
        self.rows = [
            {"filename": "late.mzML", "sample_id": "l", "sample_type": "sample"},
            {"filename": "early.mzML", "sample_id": "e", "sample_type": "QC"},
        ]
        self.scans = [
            {"file_path": str(root / "mzml/late.mzML"), "file_name": "late.mzML", "parse_status": "success",
             "acquisition_date": "2026-09-23T18:00:00+02:00", "polarity": "MS:1000129|negative scan",
             "raw_file_sha1": "a" * 40, "converted_file_sha1": "b" * 40, "checksum": "c" * 64},
            {"file_path": str(root / "mzml/early.mzML"), "file_name": "early.mzML", "parse_status": "success",
             "acquisition_date": "2026-09-23T15:00:00Z", "polarity": "MS:1000130|positive scan",
             "raw_file_sha1": "d" * 40, "converted_file_sha1": "e" * 40, "checksum": "f" * 64},
        ]
        self.save()

    def save(self):
        write_table(self.metadata, META_COLUMNS, self.rows)
        write_table(self.watcher, SCAN_COLUMNS, self.scans)

    def test_enrich_and_refresh_preserve_row_order_and_existing_values(self):
        columns, rows = helper.enrich(self.metadata, self.watcher)
        helper.save_tsv(self.metadata, columns, rows, in_place=True)
        with self.metadata.open(newline="", encoding="utf-8") as source:
            reader = csv.DictReader(source, delimiter="\t")
            self.assertEqual(reader.fieldnames[:4], ["injection_order", *META_COLUMNS])
            rows = list(reader)
        self.assertEqual([row["filename"] for row in rows], ["late.mzML", "early.mzML"])
        self.assertEqual([row["injection_order"] for row in rows], ["2", "1"])
        self.assertEqual([row["sample_id"] for row in rows], ["l", "e"])
        self.assertEqual(rows[0]["polarity"], "MS:1000129|negative scan")
        self.assertEqual(rows[0]["source_raw_sha1"], "a" * 40)
        self.assertEqual(rows[0]["converted_mzml_sha1"], "b" * 40)
        self.assertEqual(rows[0]["converted_mzml_sha256"], "c" * 64)
        self.assertEqual(rows[1]["converted_mzml_sha256"], "f" * 64)
        enriched = self.metadata.read_bytes()
        columns, rows = helper.enrich(self.metadata, self.watcher)
        helper.save_tsv(self.metadata, columns, rows, in_place=True)
        self.assertEqual(self.metadata.read_bytes(), enriched)

    def test_refresh_moves_existing_injection_order_to_front(self):
        write_table(self.metadata, [*META_COLUMNS, "injection_order"],
                    [{**row, "injection_order": "0"} for row in self.rows])
        columns, rows = helper.enrich(self.metadata, self.watcher)
        self.assertEqual(columns[:4], ["injection_order", *META_COLUMNS])
        helper.save_tsv(self.metadata, columns, rows, in_place=True)
        with self.metadata.open(newline="", encoding="utf-8") as source:
            self.assertEqual([row["injection_order"] for row in csv.DictReader(source, delimiter="\t")],
                             ["2", "1"])

    def test_unmatched_record_fails_without_changing_input(self):
        self.scans.pop()
        self.save()
        before = self.metadata.read_bytes()
        with self.assertRaisesRegex(ValueError, "No mzmlwatcher record"):
            helper.enrich(self.metadata, self.watcher)
        self.assertEqual(self.metadata.read_bytes(), before)

    def test_duplicate_watcher_filename_rejected(self):
        self.scans.append(dict(self.scans[0]))
        self.save()
        with self.assertRaisesRegex(ValueError, "Duplicate mzmlwatcher file_name"):
            helper.enrich(self.metadata, self.watcher)

    def test_tied_acquisitions_cannot_be_arbitrarily_ranked(self):
        self.scans[1]["acquisition_date"] = "2026-09-23T16:00:00Z"
        self.save()
        with self.assertRaisesRegex(ValueError, "share acquisition_date"):
            helper.enrich(self.metadata, self.watcher)

    def test_missing_input_sha_and_failed_parse_rejected(self):
        for field, value in (("raw_file_sha1", ""), ("parse_status", "error")):
            with self.subTest(field=field):
                self.scans[0][field] = value
                self.save()
                with self.assertRaisesRegex(ValueError, field):
                    helper.enrich(self.metadata, self.watcher)
                self.scans[0][field] = "a" * 40 if field == "raw_file_sha1" else "success"

    def test_missing_converted_sha256_rejected(self):
        self.scans[0]["checksum"] = ""
        self.save()
        with self.assertRaisesRegex(ValueError, "missing or invalid SHA-256 checksum"):
            helper.enrich(self.metadata, self.watcher)

    def test_watcher_record_outside_scan_directory_rejected(self):
        self.save()
        with self.assertRaisesRegex(ValueError, "outside mzML directory"):
            helper.enrich(self.metadata, self.watcher, self.metadata.parent / "wrong")

    def test_mounted_default_and_explicit_override(self):
        batch = self.metadata.parent / "group/docs/mapp_project_00016/mapp_batch_00280"
        batch.mkdir(parents=True)
        mount = self.metadata.parent / "msdata"
        expected = mount / "mapp_project_00016/mapp_batch_00280"
        expected.mkdir(parents=True)
        other = self.metadata.parent / "other_mzml"
        other.mkdir()
        with mock.patch.object(helper, "DEFAULT_MZML_ROOT", mount):
            with mock.patch.dict(os.environ, {"MZMLWATCHER_SCAN_DIR": ""}):
                self.assertEqual(helper.mzml_directory(batch, None), expected)
                self.assertEqual(helper.mzml_directory(batch, other), other)
            with mock.patch.dict(os.environ, {"MZMLWATCHER_SCAN_DIR": str(other)}):
                self.assertEqual(helper.mzml_directory(batch, None), other)
                self.assertEqual(helper.mzml_directory(batch, expected), expected)

    def test_no_mzml_directory_does_not_change_metadata(self):
        batch = self.metadata.parent / "mapp_batch_00003"
        treated = batch / "metadata/treated"
        treated.mkdir(parents=True)
        metadata = treated / f"{batch.name}_metadata.tsv"
        write_table(metadata, META_COLUMNS, self.rows)
        original = metadata.read_bytes()
        environment = os.environ.copy()
        environment.pop("MZMLWATCHER_SCAN_DIR", None)
        result = subprocess.run([sys.executable, str(SCRIPT), "--batch", str(batch)],
                                text=True, capture_output=True, env=environment)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("mzML directory not found at", result.stderr)
        self.assertEqual(metadata.read_bytes(), original)

    @unittest.skipUnless(shutil.which("mzmlwatcher"), "mzmlwatcher executable required")
    def test_batch_scans_exports_and_enriches_without_preexisting_tsv(self):
        for layout in ("in-batch", "external"):
            with self.subTest(layout=layout):
                batch = self.metadata.parent / f"mapp_batch_{'00001' if layout == 'in-batch' else '00002'}"
                treated = batch / "metadata/treated"
                treated.mkdir(parents=True)
                metadata = treated / f"{batch.name}_metadata.tsv"
                write_table(metadata, META_COLUMNS, [
                    {"filename": FIXTURE.name, "sample_id": "sample-1", "sample_type": "sample"},
                ])
                converted = batch / "mzml" if layout == "in-batch" else self.metadata.parent / "converted"
                converted.mkdir()
                mzml = converted / FIXTURE.name
                shutil.copyfile(FIXTURE, mzml)
                old = mzml.stat().st_mtime - 60
                os.utime(mzml, (old, old))
                command = [sys.executable, str(SCRIPT), "--batch", str(batch)]
                if layout == "external":
                    command.extend(["--mzml-dir", str(converted)])
                first = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
                artifacts = batch / "results/mzmlwatcher"
                self.assertTrue((artifacts / "mzmlwatcher.sqlite").is_file())
                self.assertTrue((artifacts / "mzmlwatcher.tsv").is_file())
                with metadata.open(newline="", encoding="utf-8") as source:
                    reader = csv.DictReader(source, delimiter="\t")
                    self.assertEqual(reader.fieldnames[0], "injection_order")
                    row = next(reader)
                self.assertEqual(row["injection_order"], "1")
                self.assertEqual(row["polarity"], "MS:1000130|positive scan")
                self.assertEqual(row["source_raw_sha1"], "0123456789abcdef0123456789abcdef01234567")
                self.assertEqual(row["converted_mzml_sha1"], "abcdef0123456789abcdef0123456789abcdef01")
                self.assertEqual(row["converted_mzml_sha256"], hashlib.sha256(mzml.read_bytes()).hexdigest())
                enriched = metadata.read_bytes()
                second = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
                self.assertEqual(metadata.read_bytes(), enriched)


if __name__ == "__main__":
    unittest.main()
