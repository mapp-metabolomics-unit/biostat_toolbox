#!/usr/bin/env python3

import argparse
import pandas as pd
import requests
import urllib.parse
import sys
from pathlib import Path


DEFAULT_QLEVER_ENDPOINT = "https://qlever.dev/api/wikidata"  # change if needed


def fetch_wikidata_info(inchikey: str, endpoint: str) -> dict:
    prefix = inchikey[:14]  # use only the first 14 chars (InChIKey second block)
    prefix = prefix.replace('"', '\\"')  # basic escaping for safety
    query = f"""
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT
  ?item
  (SAMPLE(?englishLabel) AS ?label)
  (GROUP_CONCAT(DISTINCT ?synonym; SEPARATOR=" | ") AS ?allSynonyms)
  (SAMPLE(?InChIValue) AS ?InChI)
  (SAMPLE(?matchingInChIKey) AS ?InChIKey)
  (SAMPLE(?canonicalSMILESValue) AS ?canonicalSMILES)
  (SAMPLE(?isomericSMILESValue) AS ?isomericSMILES)
  (SAMPLE(?massValue) AS ?mass)
WHERE {{
  ?item wdt:P235 ?matchingInChIKey .
  FILTER(STRSTARTS(?matchingInChIKey, "{prefix}"))

  OPTIONAL {{ ?item wdt:P234 ?InChIValue . }}
  OPTIONAL {{ ?item wdt:P233 ?canonicalSMILESValue . }}
  OPTIONAL {{ ?item wdt:P2017 ?isomericSMILESValue . }}
  OPTIONAL {{ ?item wdt:P2067 ?massValue . }}

  OPTIONAL {{
    ?item rdfs:label ?englishLabel .
    FILTER (LANG(?englishLabel) = "en")
  }}

  OPTIONAL {{ ?item skos:altLabel ?synonym . }}
}}
GROUP BY ?item
"""

    url = endpoint.rstrip("/") + "/api/sparql?query=" + urllib.parse.quote(query)
    response = requests.get(url)

    try:
        data = response.json()
    except Exception:
        return {}

    bindings = data.get("results", {}).get("bindings", [])
    if not bindings:
        return {}

    row = bindings[0]

    def get(field):
        return row[field]["value"] if field in row else None

    return {
        "qid": get("item"),
        "label": get("label"),
        "synonyms": get("allSynonyms"),
        "inchi": get("InChI"),
        "inchikey_returned": get("InChIKey"),
        "smiles_canonical": get("canonicalSMILES"),
        "smiles_isomeric": get("isomericSMILES"),
        "mass": get("mass"),
    }


def main():
    parser = argparse.ArgumentParser(description="Enrich a CSV of InChIKeys with Wikidata chemical metadata (via QLever SPARQL).")

    parser.add_argument("--input", "-i", required=True, help="Input CSV file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file")
    parser.add_argument("--column", "-c", required=True, help="Column name containing InChIKeys")
    parser.add_argument("--endpoint", "-e", default=DEFAULT_QLEVER_ENDPOINT,
                        help="QLever SPARQL endpoint URL (default: EarthMetabolome instance)")
    parser.add_argument("--delimiter", "-d", default=None,
                        help="Optional delimiter override (default: auto-detect based on file extension)")

    args = parser.parse_args()

    read_csv_kwargs = {}
    if args.delimiter:
        read_csv_kwargs["sep"] = args.delimiter
    else:
        suffix = Path(args.input).suffix.lower()
        if suffix == ".tsv":
            read_csv_kwargs["sep"] = "\t"
        elif suffix not in (".csv", ""):
            # Unknown extension: let pandas sniff the delimiter via python engine.
            read_csv_kwargs["sep"] = None
            read_csv_kwargs["engine"] = "python"

    try:
        df = pd.read_csv(args.input, **read_csv_kwargs)
    except Exception as exc:
        sys.exit(f"❌ Could not read input CSV: {exc}")

    if args.column not in df.columns:
        sys.exit(f"❌ Column '{args.column}' not found in CSV. Columns available: {list(df.columns)}")

    print(f"🔍 Enriching {len(df)} rows using InChIKey column '{args.column}'…")

    cache = {}  # avoid repeating identical SPARQL lookups
    results = []
    for idx, raw_value in enumerate(df[args.column]):
        print(f"  [{idx+1}/{len(df)}] {raw_value}")

        if pd.isna(raw_value):
            results.append({})
            continue

        inchikey = str(raw_value).strip()
        if not inchikey:
            results.append({})
            continue

        search_key = inchikey.upper()[:14]
        if not search_key:
            results.append({})
            continue

        if search_key not in cache:
            cache[search_key] = fetch_wikidata_info(search_key, args.endpoint)

        results.append(cache[search_key])

    enriched_df = df.join(pd.DataFrame(results))
    enriched_df.to_csv(args.output, index=False, sep="\t")

    print(f"\n✅ Done! Output written to: {args.output}")


if __name__ == "__main__":
    main()
