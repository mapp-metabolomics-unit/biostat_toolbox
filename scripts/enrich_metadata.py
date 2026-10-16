#!/usr/bin/env python3
"""Scan converted mzML files with mzmlwatcher and enrich batch metadata."""

import argparse
import csv
from datetime import datetime
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import tempfile


ADDED = (
    "injection_order",
    "acquisition_date",
    "polarity",
    "source_raw_sha1",
    "converted_mzml_sha1",
    "converted_mzml_sha256",
)
WATCHER_FIELDS = (
    "file_path", "file_name", "parse_status", "acquisition_date", "polarity",
    "raw_file_sha1", "converted_file_sha1", "checksum",
)
DEFAULT_MZML_ROOT = Path("/media/nvme8tb/msdata")
SHA1 = re.compile(r"[0-9a-fA-F]{40}\Z")
SHA256 = re.compile(r"[0-9a-fA-F]{64}\Z")


def load_tsv(path, required):
    with path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        columns = reader.fieldnames
        if not columns or len(columns) != len(set(columns)) or any(not name for name in columns):
            raise ValueError(f"{path}: missing or duplicate TSV header")
        missing = set(required) - set(columns)
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
        if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"{path}: empty or malformed TSV rows")
    return columns, rows


def enrich(metadata, watcher, scan_dir=None):
    columns, rows = load_tsv(metadata, ("filename",))
    _, scans = load_tsv(watcher, WATCHER_FIELDS)
    names = [row["filename"] for row in rows]
    wanted = set(names)
    if any(not name for name in names) or len(names) != len(wanted):
        raise ValueError("Metadata filenames must be nonempty and unique")
    selected = {}
    for scan in scans:
        name = scan["file_name"]
        if name in wanted and scan_dir is not None and not Path(scan["file_path"]).is_relative_to(scan_dir):
            raise ValueError(f"{name}: watcher record is outside mzML directory {scan_dir}")
        if name in wanted:
            if name in selected:
                raise ValueError(f"Duplicate mzmlwatcher file_name: {name}")
            selected[name] = scan
    missing = wanted - set(selected)
    if missing:
        raise ValueError(f"No mzmlwatcher record for {len(missing)} metadata filename(s): {', '.join(sorted(missing)[:5])}")

    ordered = []
    for name, scan in selected.items():
        if scan["parse_status"] != "success":
            raise ValueError(f"{name}: mzmlwatcher parse_status is {scan['parse_status']!r}")
        try:
            date = datetime.fromisoformat(scan["acquisition_date"].replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError(f"{name}: invalid acquisition_date") from exc
        if date.tzinfo is None or date.utcoffset() is None:
            raise ValueError(f"{name}: acquisition_date must include a timezone")
        for field in ("raw_file_sha1", "converted_file_sha1"):
            if not SHA1.fullmatch(scan[field]):
                raise ValueError(f"{name}: missing or invalid {field}")
        if not SHA256.fullmatch(scan["checksum"]):
            raise ValueError(f"{name}: missing or invalid SHA-256 checksum")
        if not scan["polarity"]:
            raise ValueError(f"{name}: missing polarity")
        ordered.append((date, name))
    ordered.sort()
    for previous, current in zip(ordered, ordered[1:]):
        if previous[0] == current[0]:
            raise ValueError(f"Cannot determine injection order: {previous[1]} and {current[1]} share acquisition_date")
    ranks = {name: str(index) for index, (_, name) in enumerate(ordered, 1)}
    for row in rows:
        scan = selected[row["filename"]]
        row.update({
            "injection_order": ranks[row["filename"]],
            "acquisition_date": scan["acquisition_date"],
            "polarity": scan["polarity"],
            "source_raw_sha1": scan["raw_file_sha1"],
            "converted_mzml_sha1": scan["converted_file_sha1"],
            "converted_mzml_sha256": scan["checksum"],
        })
    remaining = columns + [field for field in ADDED if field not in columns]
    return ["injection_order"] + [field for field in remaining if field != "injection_order"], rows


def save_tsv(path, columns, rows, *, in_place):
    temp = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=path.parent,
                                         prefix=f".{path.name}.", delete=False) as target:
            temp = Path(target.name)
            writer = csv.DictWriter(target, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        if in_place:
            os.chmod(temp, stat.S_IMODE(path.stat().st_mode))
            os.replace(temp, path)
        else:
            os.link(temp, path)  # Atomic publication; never overwrite an existing output.
    finally:
        if temp is not None:
            temp.unlink(missing_ok=True)


def mzml_directory(batch, specified):
    default = DEFAULT_MZML_ROOT / batch.parent.name / batch.name
    if specified is not None:
        directory = specified
    elif os.environ.get("MZMLWATCHER_SCAN_DIR"):
        directory = Path(os.environ["MZMLWATCHER_SCAN_DIR"])
    elif default.is_dir():
        directory = default
    elif (batch / "mzml").is_dir():
        directory = batch / "mzml"
    elif any(path.is_file() and path.suffix.lower() == ".mzml" for path in batch.iterdir()):
        directory = batch
    else:
        raise ValueError(f"mzML directory not found at {default}; pass --mzml-dir or set MZMLWATCHER_SCAN_DIR")
    directory = directory.resolve()
    if not directory.is_dir():
        raise ValueError(f"mzML directory does not exist: {directory}")
    return directory


def scan_mzml(batch, directory):
    artifacts = batch / "results" / "mzmlwatcher"
    artifacts.mkdir(parents=True, exist_ok=True)
    sqlite = artifacts / "mzmlwatcher.sqlite"
    tsv = artifacts / "mzmlwatcher.tsv"
    subprocess.run(["mzmlwatcher", "scan", str(directory), "--sqlite", str(sqlite),
                    "--tsv", str(tsv), "--recursive", "--checksum", "sha256"], check=True)
    return tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", required=True, type=Path, help="batch directory with treated metadata")
    parser.add_argument("--mzml-dir", type=Path, help="override /media/nvme8tb/msdata/<project>/<batch>")
    parser.add_argument("--output", type=Path, help="write a new metadata TSV instead of updating in place")
    args = parser.parse_args()
    try:
        batch = args.batch.resolve()
        if not batch.is_dir():
            raise ValueError(f"Batch directory does not exist: {batch}")
        metadata = batch / "metadata" / "treated" / f"{batch.name}_metadata.tsv"
        if metadata.is_symlink():
            raise ValueError("Refusing to replace a metadata symlink")
        load_tsv(metadata, ("filename",))  # Fail before a potentially expensive scan.
        if args.output and args.output.resolve() == metadata:
            raise ValueError("Omit --output to update metadata in place")
        if args.output and args.output.exists():
            raise ValueError(f"Output already exists: {args.output}")
        directory = mzml_directory(batch, args.mzml_dir)
        watcher_tsv = scan_mzml(batch, directory)
        columns, rows = enrich(metadata, watcher_tsv, directory)
        destination = args.output or metadata
        save_tsv(destination, columns, rows, in_place=args.output is None)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"enrich_metadata: {exc}\n")
    print(f"Enriched {len(rows)} injections: {destination}")


if __name__ == "__main__":
    main()
