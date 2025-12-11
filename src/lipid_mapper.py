"""
Utility for mapping SIRIUS structure identifications to SwissLipids IDs or
LipidMaps metadata.

The tool expects a TSV file such as ``structure_identifications.tsv`` created by
SIRIUS.  For every row it extracts ``InChIkey2D`` and ``mappingFeatureId``, runs
the selected lipid database lookup and writes a compact TSV with the feature id,
the 2D InChIKey and the retrieved identifiers/metadata.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

import requests
from tqdm import tqdm

SWISSLIPIDS_SEARCH_URL = "https://www.swisslipids.org/api/index.php/advancedSearch"
LIPIDMAPS_URL_TEMPLATE = (
    "https://www.lipidmaps.org/rest/compound/inchi_key/{inchikey}/all"
)
BASE_HEADERS = ("mappingFeatureId", "InChIkey2D")
SWISS_HEADERS = ("SwissLipidsID",)
LIPIDMAPS_FIELDS = [
    "input",
    "regno",
    "lm_id",
    "name",
    "sys_name",
    "synonyms",
    "abbrev",
    "core",
    "main_class",
    "sub_class",
    "exactmass",
    "formula",
    "inchi",
    "inchi_key",
    "kegg_id",
    "hmdb_id",
    "chebi_id",
    "lipidbank_id",
    "pubchem_cid",
    "smiles",
]
LIPIDMAPS_HEADERS = tuple(f"lipid_maps_{field}" for field in LIPIDMAPS_FIELDS)


def normalize_field_lookup(row: Dict[str, str], target: str) -> Optional[str]:
    """Return a field value ignoring column name casing, if present."""
    target_lower = target.lower()
    for column, value in row.items():
        if column and column.lower() == target_lower:
            return value
    return None


def log_warning(message: str, *, progress=None) -> None:
    if progress:
        progress.write(f"[warning] {message}")
    else:
        print(f"[warning] {message}", file=sys.stderr)


def log_info(message: str, *, progress=None) -> None:
    if progress:
        progress.write(f"[info] {message}")
    else:
        print(f"[info] {message}", file=sys.stderr)


def fetch_swisslipids_data(
    inchikey: str, *, session: requests.Session, retries: int = 3, progress=None
) -> Optional[Dict[str, str]]:
    """Run the SwissLipids advanced search API for a given InChIKey."""
    params = {"InChIkey": inchikey}
    for attempt in range(1, retries + 1):
        try:
            response = session.get(SWISSLIPIDS_SEARCH_URL, params=params, timeout=30)
            response.raise_for_status()
            payload = response.json()
        except (requests.RequestException, ValueError) as exc:
            if attempt >= retries:
                log_warning(
                    f"SwissLipids lookup failed for {inchikey}: {exc}", progress=progress
                )
                return None
            continue

        if isinstance(payload, list) and payload:
            entity = payload[0]
            swiss_id = entity.get("entity_id")
            if swiss_id:
                return {"SwissLipidsID": swiss_id}
        return None
    return None


def fetch_lipidmaps_data(
    inchikey: str, *, session: requests.Session, progress=None
) -> Optional[Dict[str, str]]:
    """Query LipidMaps REST API for a given InChIKey and prefix the fields."""
    inchikey_full = ensure_standard_inchikey(inchikey)
    url = LIPIDMAPS_URL_TEMPLATE.format(inchikey=inchikey_full)
    try:
        response = session.get(url, timeout=30)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        log_warning(f"LipidMaps lookup failed for {inchikey}: {exc}", progress=progress)
        return None

    if not isinstance(payload, dict):
        return None

    normalized = {}
    for field in LIPIDMAPS_FIELDS:
        value = payload.get(field)
        key = f"lipid_maps_{field}"
        normalized[key] = "" if value is None else str(value)
    return normalized


def output_headers_for_source(source: str) -> Iterable[str]:
    if source == "swisslipids":
        return BASE_HEADERS + SWISS_HEADERS
    if source == "lipidmaps":
        return BASE_HEADERS + LIPIDMAPS_HEADERS
    raise ValueError(f"Unsupported source: {source}")


def fetcher_for_source(source: str):
    if source == "swisslipids":
        return fetch_swisslipids_data
    if source == "lipidmaps":
        return fetch_lipidmaps_data
    raise ValueError(f"Unsupported source: {source}")


def map_structures(input_path: Path, output_path: Path, source: str) -> None:
    """Read the TSV, query the selected API, and produce the mapping TSV."""
    cache: Dict[str, Optional[Dict[str, str]]] = {}
    total_rows = 0
    mapped = 0

    total_expected = count_data_rows(input_path)
    output_headers = output_headers_for_source(source)
    fetch_entry = fetcher_for_source(source)

    with input_path.open("r", encoding="utf-8", newline="") as handle_in, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as handle_out:
        reader = csv.DictReader(handle_in, delimiter="\t")
        writer = csv.DictWriter(handle_out, delimiter="\t", fieldnames=output_headers)
        writer.writeheader()

        session = requests.Session()
        progress_label = f"{source.capitalize()} lookups"
        progress = tqdm(
            total=total_expected or None,
            desc=progress_label,
            unit="row",
        )
        for row in reader:
            total_rows += 1
            inchikey = normalize_field_lookup(row, "InChIkey2D")
            feature_id = normalize_field_lookup(row, "mappingFeatureId")

            if not inchikey or not feature_id:
                log_warning(
                    f"skipped row {total_rows}: missing InChIkey2D or mappingFeatureId",
                    progress=progress,
                )
                progress.update(1)
                continue

            if inchikey not in cache:
                cache[inchikey] = fetch_entry(
                    inchikey, session=session, progress=progress
                )
            metadata = cache[inchikey]
            if not metadata:
                progress.update(1)
                continue
            mapped += 1

            base_entry = {"mappingFeatureId": feature_id, "InChIkey2D": inchikey}
            base_entry.update(metadata)
            writer.writerow(base_entry)
            progress.update(1)

        progress.close()

    log_info(
        f"processed {total_rows} rows -> {mapped} {source} hits",
        progress=None,
    )


def count_data_rows(path: Path) -> int:
    """Return the number of non-header lines in a TSV file."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        total = sum(1 for _ in handle)
    return max(total - 1, 0)


def ensure_standard_inchikey(inchikey: str) -> str:
    """Append the standard key suffix if missing (required by LipidMaps endpoint)."""
    suffix = "-UHFFFAOYSA-N"
    if inchikey.endswith(suffix):
        return inchikey
    base = inchikey.split("-")[0]
    return f"{base}{suffix}"


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map SIRIUS structure_identifications.tsv entries to lipid databases",
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to the SIRIUS structure_identifications.tsv file",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path for the output TSV (mappingFeatureId, InChIkey2D, plus API fields)",
    )
    parser.add_argument(
        "--source",
        choices=("swisslipids", "lipidmaps"),
        default="swisslipids",
        help="Which API to query for mappings",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Iterable[str]] = None) -> None:
    args = parse_args(argv)
    if not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")
    map_structures(args.input, args.output, args.source)


if __name__ == "__main__":
    main()
