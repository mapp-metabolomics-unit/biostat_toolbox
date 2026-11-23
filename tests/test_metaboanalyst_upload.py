"""Consumer-visible boundaries of the MetaboAnalyst web exporter."""

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


EXPORTER = Path(__file__).resolve().parents[1] / "scripts" / "prepare_metaboanalyst_upload.py"
spec = importlib.util.spec_from_file_location("prepare_metaboanalyst_upload", EXPORTER)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def write_tsv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.reader(source))


class MetaboAnalystUploadTests(unittest.TestCase):
    def test_upload_excludes_qc_and_incomplete_pairs_without_imputing_peaks(self):
        with tempfile.TemporaryDirectory() as tmp:
            study = Path(tmp) / "study"
            tables = study / "tables"
            tables.mkdir(parents=True)
            (study / "COMPLETE").write_text("complete", encoding="utf-8")
            (study / "manifest.yaml").write_text("kind: study\n", encoding="utf-8")
            header = ["filename", "sample_id", "sample_type", "ATTRIBUTE_plant_uid",
                      "ATTRIBUTE_group", "ATTRIBUTE_timepoint", "ATTRIBUTE_batch"]
            samples, intensities = [], []
            for group in ("Infected+treated", "Infected_control"):
                for plant in range(3):
                    for time in (0, 1):
                        sample = f"S_{group.replace('+', '_')}_{plant}_T{time}"
                        samples.append([f"{sample}.mzML", sample, "sample", f"{group}_{plant}", group, f"T{time}", "B1"])
                        intensities.append([f"{sample}.mzML", "0" if plant == 0 and time == 0 else "12.5", "8", "0"])
            samples.append(["S_incomplete.mzML", "S_incomplete", "sample", "incomplete", "Infected_control", "T0", "B1"])
            intensities.append(["S_incomplete.mzML", "99", "14", "7"])
            samples.extend([["QC.mzML", "QC", "QC", "", "", "", ""], ["blank.mzML", "blank", "blank", "", "", "", ""]])
            intensities.extend([["QC.mzML", "32", "8", "0"], ["blank.mzML", "200", "3", "0"]])
            write_tsv(tables / "sample_metadata.tsv", header, samples)
            write_tsv(tables / "variable_metadata.tsv", ["feature_id", "row m/z", "row retention time"],
                      [["4", "227.9", "0.47"], ["6", "249.9", "0.48"], ["8", "251.5", "0.51"]])
            write_tsv(tables / "matrix_raw.tsv", ["sample", "4", "6", "8"], intensities)
            output = Path(tmp) / "upload"
            module.export(study, output)
            t0 = read_csv(output / "statistical_T0.csv")
            t1 = read_csv(output / "statistical_T1.csv")
            metadata = read_csv(output / "time_series_metadata.csv")
            series = read_csv(output / "time_series.csv")
            self.assertEqual([len(t0), len(t1), len(series), len(metadata)], [8, 7, 13, 13])
            self.assertEqual(t0[0], ["Sample", "Class", "F_4", "F_6", "F_8"])
            self.assertEqual(t1[0], ["Sample", "Class", "F_4", "F_6"])
            self.assertEqual(series[0], ["Sample", "F_4", "F_6"])
            self.assertEqual(metadata[0], ["Sample", "Phenotype", "Time", "Subject"])
            self.assertEqual(next(row for row in t0 if row[0] == "S_Infected_treated_0_T0")[2], "")
            self.assertNotIn("S_incomplete", [row[0] for row in series[1:]])
            self.assertEqual([row[0] for row in series[1:]], [row[0] for row in metadata[1:]])
            self.assertFalse(any(row[0] in ("QC", "blank") for row in t0[1:] + t1[1:] + series[1:]))


if __name__ == "__main__":
    unittest.main()
