use std::collections::HashSet;
use std::fs::{self, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};
use serde_json::Value;

#[derive(Clone, PartialEq, Eq)]
pub struct Table {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

fn reader_without_bom(path: &Path) -> Result<impl Read> {
    let file = fs::File::open(path).with_context(|| format!("Opening {}", path.display()))?;
    let mut reader = BufReader::new(file);
    // Strip only an initial UTF-8 BOM, leaving every actual cell unchanged.
    if reader.fill_buf()?.starts_with(&[0xef, 0xbb, 0xbf]) {
        reader.consume(3);
    }
    Ok(reader)
}

fn checked_header(
    header: &csv::StringRecord,
    label: &str,
    allow_empty: bool,
) -> Result<Vec<String>> {
    let columns: Vec<String> = header.iter().map(str::to_owned).collect();
    let mut seen = HashSet::new();
    // MZmine may emit a single trailing unnamed column; it is not an intensity field.
    if columns.is_empty()
        || columns
            .iter()
            .any(|column| (!allow_empty && column.trim().is_empty()) || !seen.insert(column))
    {
        bail!(
            "{label} header must have unique{} columns",
            if allow_empty { "" } else { ", nonempty" }
        );
    }
    Ok(columns)
}

impl Table {
    pub fn read_tsv(path: &Path) -> Result<Self> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(reader_without_bom(path)?);
        let columns = checked_header(reader.headers()?, "Metadata", false)?;
        let mut rows = Vec::new();
        for result in reader.records() {
            let record = result.with_context(|| format!("Reading {}", path.display()))?;
            if record.len() != columns.len() || record.iter().all(str::is_empty) {
                bail!(
                    "Metadata must have nonempty rows with {} columns",
                    columns.len()
                );
            }
            rows.push(record.iter().map(str::to_owned).collect());
        }
        if rows.is_empty() {
            bail!("Metadata must have at least one row");
        }
        Ok(Self { columns, rows })
    }

    pub fn column(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|column| column == name)
    }

    pub fn write_tsv(&self, path: &Path) -> Result<()> {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .with_context(|| format!("Creating {}", path.display()))?;
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .terminator(csv::Terminator::Any(b'\n'))
            .from_writer(file);
        writer.write_record(&self.columns)?;
        for row in &self.rows {
            if row.len() != self.columns.len() {
                bail!("Metadata row does not match its header width");
            }
            writer.write_record(row)?;
        }
        writer.flush()?;
        Ok(())
    }
}

pub fn quant_header(path: &Path) -> Result<Vec<String>> {
    let mut reader = csv::ReaderBuilder::new().from_reader(reader_without_bom(path)?);
    let columns = checked_header(reader.headers()?, "Quantification", true)?;
    if !columns.iter().any(|column| column == "row ID") {
        bail!("Quantification lacks the required feature column \"row ID\"");
    }
    Ok(columns)
}

fn write_yaml(path: &Path, value: &Value) -> Result<()> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .with_context(|| format!("Creating {}", path.display()))?;
    // JSON is valid YAML; quoting every string keeps YAML 1.1 readers from
    // treating metadata labels such as "yes" and "no" as booleans.
    serde_json::to_writer_pretty(&mut file, value)?;
    file.write_all(b"\n")?;
    file.flush()?;
    Ok(())
}

pub fn save(
    folder: &Path,
    dataset: &Value,
    recipe: &Value,
    metadata: Option<&Table>,
) -> Result<()> {
    let parent = folder
        .parent()
        .context("Recipe folder needs a parent directory")?;
    fs::create_dir_all(parent).with_context(|| format!("Creating {}", parent.display()))?;
    // create_dir is the exclusive claim: an existing recipe is never touched, even on failure.
    fs::create_dir(folder).with_context(|| format!("Creating {}", folder.display()))?;
    let result = (|| {
        if let Some(table) = metadata {
            table.write_tsv(&folder.join("metadata.tsv"))?;
        }
        write_yaml(&folder.join("dataset.yaml"), dataset)?;
        write_yaml(&folder.join("recipe.yaml"), recipe)
    })();
    if let Err(error) = result {
        fs::remove_dir_all(folder).with_context(|| {
            format!(
                "Removing incomplete recipe {} after: {error}",
                folder.display()
            )
        })?;
        return Err(error);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mzmine_trailing_unnamed_column_is_not_a_missing_sample() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let csv = dir.path().join("quant.csv");
        fs::write(
            &csv,
            "row ID,A.mzML Peak height,B.mzML Peak height,\n1,10,20,\n",
        )?;
        assert_eq!(
            quant_header(&csv)?,
            ["row ID", "A.mzML Peak height", "B.mzML Peak height", ""]
        );
        fs::write(
            &csv,
            "row ID,A.mzML Peak height,A.mzML Peak height,\n1,10,20,\n",
        )?;
        assert!(quant_header(&csv).is_err());
        Ok(())
    }

    #[test]
    fn existing_recipe_directory_is_never_overwritten() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let folder = dir.path().join("recipes").join("new-recipe");
        let metadata = Table {
            columns: vec!["filename".into(), "sample_type".into()],
            rows: vec![vec!["blank_new.mzML".into(), "blank".into()]],
        };
        let dataset = serde_json::json!({"metadata": "recipes/new-recipe/metadata.tsv"});
        save(
            &folder,
            &dataset,
            &serde_json::json!({"seed": 1}),
            Some(&metadata),
        )?;
        assert!(fs::read_to_string(folder.join("metadata.tsv"))?.contains("blank_new.mzML"));
        let original = fs::read(folder.join("recipe.yaml"))?;
        assert!(save(
            &folder,
            &dataset,
            &serde_json::json!({"seed": 2}),
            Some(&metadata)
        )
        .is_err());
        assert_eq!(fs::read(folder.join("recipe.yaml"))?, original);
        Ok(())
    }

    #[test]
    fn saved_configs_keep_metadata_labels_as_strings_for_yaml_readers() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let folder = dir.path().join("recipes").join("treatment");
        let dataset = serde_json::json!({
            "metadata": "2026-01-01.tsv",
            "roles": {"levels": {"sample": "on", "qc": "off", "blank": "null"}}
        });
        let recipe = serde_json::json!({
            "sample_plan": {"groups": [
                {"source": "no", "target": "no"},
                {"source": "yes", "target": "yes"}
            ]},
            "analyses": {"pca": {"enabled": true}},
            "seed": 1
        });
        save(&folder, &dataset, &recipe, None)?;
        let saved_dataset: Value =
            serde_json::from_reader(fs::File::open(folder.join("dataset.yaml"))?)?;
        let saved_recipe: Value =
            serde_json::from_reader(fs::File::open(folder.join("recipe.yaml"))?)?;
        assert_eq!(saved_dataset, dataset);
        assert_eq!(saved_recipe, recipe);
        Ok(())
    }
}
