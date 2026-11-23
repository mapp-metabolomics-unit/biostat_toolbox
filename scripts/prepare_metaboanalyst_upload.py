#!/usr/bin/env python3
"""Prepare MetaboAnalyst 6.0 web uploads from one sealed biostat_toolbox study."""

import argparse
from collections import Counter
import csv
import math
from pathlib import Path
import re


SAFE_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")
TIMEPOINT = re.compile(r"T(\d+)\Z")


def rows(path, required):
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing required columns {set(required) - set(reader.fieldnames or ())}")
        result = list(reader)
    if not result or any(None in row or None in row.values() for row in result):
        raise ValueError(f"{path}: empty or malformed TSV")
    return result


def clean_label(label):
    return re.sub(r"[^A-Za-z0-9_]", "_", label)


def export(study, output):
    tables = study / "tables"
    if not (study / "COMPLETE").is_file() or not (study / "manifest.yaml").is_file():
        raise ValueError(f"{study}: not a completed sealed study")
    metadata = rows(tables / "sample_metadata.tsv", (
        "filename", "sample_id", "sample_type", "ATTRIBUTE_plant_uid", "ATTRIBUTE_group",
        "ATTRIBUTE_timepoint", "ATTRIBUTE_batch",
    ))
    features = rows(tables / "variable_metadata.tsv", ("feature_id", "row m/z", "row retention time"))
    feature_ids = [row["feature_id"] for row in features]
    aliases = ["F_" + re.sub(r"[^A-Za-z0-9_]", "_", fid) for fid in feature_ids]
    if any(not fid for fid in feature_ids) or len(set(feature_ids)) != len(features) or len(set(aliases)) != len(features):
        raise ValueError("Retained feature IDs must be present and distinct after safe naming")

    samples = [row for row in metadata if row["sample_type"] == "sample"]
    if not samples:
        raise ValueError("No biological samples in sealed study")
    filenames = [row["filename"] for row in samples]
    sample_ids = [row["sample_id"] for row in samples]
    if (len(set(filenames)) != len(samples) or len(set(sample_ids)) != len(samples)
            or any(not SAFE_NAME.fullmatch(sid) for sid in sample_ids)):
        raise ValueError("Biological filenames and upload sample IDs must be unique and web-safe")
    groups = {row["ATTRIBUTE_group"] for row in samples}
    group_aliases = {group: clean_label(group) for group in groups}
    if (not groups or any(not SAFE_NAME.fullmatch(alias) for alias in group_aliases.values())
            or len(set(group_aliases.values())) != len(groups)):
        raise ValueError("Biological groups cannot be uniquely represented by web-safe labels")
    for row in samples:
        if (not row["ATTRIBUTE_plant_uid"] or not row["ATTRIBUTE_batch"]
                or not TIMEPOINT.fullmatch(row["ATTRIBUTE_timepoint"])):
            raise ValueError(f"{row['sample_id']}: missing subject, batch, or T-number timepoint")
    by_filename = {row["filename"]: row for row in samples}
    matrix = {}
    with (tables / "matrix_raw.tsv").open(encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        if header[0] != "sample" or len(header) != len(set(header)) or not set(feature_ids).issubset(header):
            raise ValueError("Raw matrix has missing or duplicate feature IDs")
        positions = [header.index(fid) for fid in feature_ids]
        for line in reader:
            if len(line) != len(header):
                raise ValueError("Malformed raw-matrix row")
            filename = line[0]
            if filename not in by_filename:
                continue  # QC and blank injections are not biological groups.
            if filename in matrix:
                raise ValueError(f"Duplicate raw-matrix sample: {filename}")
            values = [line[index] for index in positions]
            for index, value in enumerate(values):
                if value in ("", "NA", "NaN"):
                    values[index] = ""
                    continue
                try:
                    number = float(value)
                except ValueError as exc:
                    raise ValueError(f"{filename}: nonnumeric intensity for {feature_ids[index]}") from exc
                if not math.isfinite(number) or number < 0:
                    raise ValueError(f"{filename}: nonfinite or negative intensity for {feature_ids[index]}")
                if number == 0:
                    values[index] = ""  # Zero from peak picking is missing, not an observed intensity.
            matrix[filename] = values
    if set(matrix) != set(filenames):
        raise ValueError(f"Raw matrix is missing {len(set(filenames) - set(matrix))} biological samples")

    by_subject = {}
    for row in samples:
        by_subject.setdefault(row["ATTRIBUTE_plant_uid"], []).append(row)
    timepoints = sorted({int(TIMEPOINT.fullmatch(row["ATTRIBUTE_timepoint"]).group(1)) for row in samples})
    complete = set()
    for subject, subject_rows in by_subject.items():
        subject_times = [int(TIMEPOINT.fullmatch(row["ATTRIBUTE_timepoint"]).group(1)) for row in subject_rows]
        if len(subject_times) != len(set(subject_times)):
            raise ValueError(f"{subject}: multiple biological injections for the same time point")
        if len({row["ATTRIBUTE_group"] for row in subject_rows}) != 1 or len({row["ATTRIBUTE_batch"] for row in subject_rows}) != 1:
            raise ValueError(f"{subject}: phenotype or batch changes during the time series")
        if len(subject_rows) == len(timepoints) and set(subject_times) == set(timepoints):
            complete.add(subject)
    if not complete or len(timepoints) < 2:
        raise ValueError("No balanced repeated-subject time series")

    def write_csv(name, header, data):
        with (output / name).open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target)
            writer.writerow(header)
            writer.writerows(data)

    output.mkdir(parents=True, exist_ok=True)
    per_time_features = {}
    for time in timepoints:
        subset = sorted((row for row in samples if row["ATTRIBUTE_timepoint"] == f"T{time}"),
                        key=lambda row: (group_aliases[row["ATTRIBUTE_group"]], row["ATTRIBUTE_plant_uid"]))
        if min(Counter(group_aliases[row["ATTRIBUTE_group"]] for row in subset).values()) < 3:
            raise ValueError(f"T{time}: fewer than three biological replicates in a group")
        available = [index for index in range(len(aliases))
                     if any(matrix[row["filename"]][index] for row in subset)]
        if not available:
            raise ValueError(f"T{time}: no observed biological features")
        per_time_features[time] = len(available)
        write_csv(f"statistical_T{time}.csv", ["Sample", "Class", *(aliases[index] for index in available)], (
            [row["sample_id"], group_aliases[row["ATTRIBUTE_group"]],
             *(matrix[row["filename"]][index] for index in available)] for row in subset
        ))
    paired = sorted((row for row in samples if row["ATTRIBUTE_plant_uid"] in complete),
                    key=lambda row: (row["ATTRIBUTE_plant_uid"], int(TIMEPOINT.fullmatch(row["ATTRIBUTE_timepoint"]).group(1))))
    paired_features = [index for index in range(len(aliases))
                       if any(matrix[row["filename"]][index] for row in paired)]
    if not paired_features:
        raise ValueError("No observed biological features in complete plants")
    write_csv("time_series.csv", ["Sample", *(aliases[index] for index in paired_features)], (
        [row["sample_id"], *(matrix[row["filename"]][index] for index in paired_features)] for row in paired
    ))
    write_csv("time_series_metadata.csv", ["Sample", "Phenotype", "Time", "Subject"], (
        [row["sample_id"], group_aliases[row["ATTRIBUTE_group"]],
         TIMEPOINT.fullmatch(row["ATTRIBUTE_timepoint"]).group(1),
         row["ATTRIBUTE_plant_uid"]] for row in paired
    ))
    write_csv("feature_map.csv", ["UploadFeature", "SourceFeature", "mz", "rt"], (
        [alias, row["feature_id"], row["row m/z"], row["row retention time"]]
        for alias, row in zip(aliases, features)
    ))
    write_csv("sample_map.csv", ["Sample", "SourceFile", "SourceGroup", "UploadClass", "Time", "Subject", "Batch", "CompleteSubject"], (
        [row["sample_id"], row["filename"], row["ATTRIBUTE_group"], group_aliases[row["ATTRIBUTE_group"]],
         row["ATTRIBUTE_timepoint"], row["ATTRIBUTE_plant_uid"], row["ATTRIBUTE_batch"],
         "yes" if row["ATTRIBUTE_plant_uid"] in complete else "no"] for row in samples
    ))
    guide = f"""# MetaboAnalyst online statistical uploads (6.0 / developer 7.0)

Source: sealed biostat_toolbox study `{study.name}`. Raw MZmine peak **heights** for {len(features)} QC/blank-filtered features; {len(samples)} biological samples, no QC or blanks. Zero/absent peaks are empty (missing); no imputation, normalization, transformation, or scaling was applied to these upload files. MetaboAnalyst's downstream choices may therefore differ from the biostat_toolbox dashboard.

Features with **no measured positive height** in a given upload cohort are omitted from that CSV only: {', '.join(f'T{time}: {per_time_features[time]}' for time in timepoints)}; matched time series: {len(paired_features)}. `feature_map.csv` covers the complete sealed feature set. Within-file sparse peaks may be removed by MetaboAnalyst's filtering; inspect its data integrity report.

- **Statistical Analysis [one factor]:** on the [developer 7.0 upload screen](https://dev.metaboanalyst.ca/MetaboAnalyst/upload/StatUploadView.xhtml), under **A plain text file (.txt or .csv)** select **Data Type: Peak intensities**, **Format: Samples in rows (unpaired)**, **Data File:** {', '.join(f'`statistical_T{time}.csv`' for time in timepoints)} **one file at a time**, then **Submit**. Each file contains one time point (independent plants); class names are in column 2. **Do not select paired:** the paired option requires a different two-group file with signed integer pair labels, which these treatment-group CSVs do not contain. Do not concatenate time points for a one-factor t-test/ANOVA: plants recur.
- **Statistical Analysis [multi factors / metadata table]:** for repeated plants, leave the single-factor upload page and use the separate metadata-table module. Upload `time_series.csv` (**samples in rows**) plus `time_series_metadata.csv`, choosing the time-series + one experimental factor design: **Phenotype** = treatment group, **Time** = {'/'.join(str(time) for time in timepoints)}, **Subject** = matched plant ID. {len(complete)} complete plants ({len(paired)} injections) are included; {len(samples) - len(paired)} incomplete-plant injections are omitted only here. This specialized 3-factor design does **not** include a batch adjustment; batch is recorded in `sample_map.csv` for inspection. Do not equate its results with batch-adjusted biostat_toolbox analyses. Keep each subject's time points together.
- `feature_map.csv` maps web-safe feature names to the original MZmine ID, m/z and retention time. `sample_map.csv` maps web-safe sample IDs and class aliases back to the source data. These maps are **not upload files**.
- Do not interpret numeric LC-MS features as confirmed compounds or upload these tables as annotated metabolite sets. Do not add QC/blank injections as biological class members.
- **Functional Analysis [MS Peaks to Pathways]** takes a different m/z-and-retention-time peak table or a contrast-specific peak list with statistical scores; the `F_` feature labels in these statistical uploads are not valid m/z–RT identifiers for that module. No pathway or biomarker training file is generated.

To regenerate from the batch directory (toolbox and batch repositories side by side):

```sh
python3 ../../../../biostat_toolbox/scripts/prepare_metaboanalyst_upload.py \\
  --study results/stats/{study.name} --output results/metaboanalyst
```

Format reference: https://www.metaboanalyst.ca/docs/Format.xhtml (Statistical Analysis [one factor] and [metadata table]). No automatic upload, online validation or consent to transfer this study's data is performed.
"""
    (output / "README.md").write_text(guide, encoding="utf-8")
    print(f"Wrote {len(timepoints)} one-factor tables and {len(paired)} matched observations for {len(features)} features to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True, help="Complete sealed study directory under results/stats")
    parser.add_argument("--output", type=Path, required=True, help="Directory for MetaboAnalyst upload files")
    arguments = parser.parse_args()
    export(arguments.study, arguments.output)
