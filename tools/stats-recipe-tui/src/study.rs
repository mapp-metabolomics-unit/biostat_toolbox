use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::path::{Component, Path, PathBuf};

use anyhow::{bail, Context, Result};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

fn hash(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}

fn regular_file(path: &Path) -> Result<()> {
    if !fs::symlink_metadata(path)?.file_type().is_file() {
        bail!("Expected a regular, non-symlink file: {}", path.display());
    }
    Ok(())
}

fn checksum(path: &Path) -> Result<String> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut sha = Sha256::new();
    let mut block = [0u8; 65536];
    loop {
        let count = reader.read(&mut block)?;
        if count == 0 { break; }
        sha.update(&block[..count]);
    }
    Ok(format!("{:x}", sha.finalize()))
}

fn safe_relative(path: &str) -> bool {
    !path.contains('\\')
        && path.split('/').all(|part| !part.is_empty() && part != "." && part != "..")
        && Path::new(path).components().all(|part| matches!(part, Component::Normal(_)))
}

fn checked_source(batch: &Path, original_batch: &Path, absolute: &str) -> Result<(PathBuf, String)> {
    let path = Path::new(absolute);
    if !path.is_absolute() || !original_batch.is_absolute() {
        bail!("Study source or batch path is not absolute: {absolute}");
    }
    let relative = path.strip_prefix(original_batch)
        .context("Study source lies outside its original batch")?
        .to_str().context("Study source path must be UTF-8")?.to_owned();
    if !safe_relative(&relative) { bail!("Invalid study source path: {relative}"); }
    let relocated = batch.join(&relative);
    let canonical = relocated.canonicalize()
        .with_context(|| format!("Missing relocated study source {}", relocated.display()))?;
    if relocated != canonical { bail!("Study source path uses a symlink or noncanonical component: {}", relocated.display()); }
    // A symlink in any path component would make the source silently mutable by redirection.
    let mut cursor = batch.to_path_buf();
    for part in Path::new(&relative).components() {
        cursor.push(part);
        if fs::symlink_metadata(&cursor)?.file_type().is_symlink() {
            bail!("Study source contains a symlink: {}", cursor.display());
        }
    }
    regular_file(&canonical)?;
    Ok((canonical, relative))
}

pub struct Study {
    pub hash: String,
    pub dataset: Value,
    pub metadata_path: PathBuf,
    pub quant_path: PathBuf,
}

impl Study {
    pub fn load(batch: &Path, study_hash: &str) -> Result<Self> {
        if !hash(study_hash) { bail!("--study must be a full lowercase 64-character SHA-256 hash"); }
        let directory = batch.join("results/stats").join(study_hash);
        for name in ["results", "results/stats"] {
            let parent = batch.join(name);
            if !fs::symlink_metadata(&parent)?.file_type().is_dir() {
                bail!("Study results path is not a real directory: {}", parent.display());
            }
        }
        if !fs::symlink_metadata(&directory)?.file_type().is_dir() {
            bail!("Study run is not a real directory: {}", directory.display());
        }
        for name in ["COMPLETE", "manifest.yaml"] { regular_file(&directory.join(name))?; }
        let marker: Value = serde_yaml::from_reader(File::open(directory.join("COMPLETE"))?)?;
        let manifest: Value = serde_yaml::from_reader(File::open(directory.join("manifest.yaml"))?)?;
        if marker["run_hash"] != study_hash || !marker["content_sha256"].as_str().is_some_and(hash)
            || manifest["kind"] != "study" || manifest["schema_version"] != "1.0.0"
            || manifest["status"] != "complete" || manifest["hashes"]["run_hash"] != study_hash
            || manifest["effective_recipe"]["kind"] != "study"
            || !manifest["effective_recipe"]["preprocessing"].is_object()
            || !manifest["effective_recipe"]["analyses"].is_object()
            || !manifest["hashes"]["dataset_hash"].as_str().is_some_and(hash)
            || !manifest["hashes"]["preprocess_hash"].as_str().is_some_and(hash) {
            bail!("Invalid sealed study marker or manifest: {}", directory.display());
        }
        let outputs = manifest["outputs"].as_object().context("Study manifest has no output map")?;
        if ["tables/sample_metadata.tsv", "objects/dataset.rds", "objects/processed.rds"]
            .iter().any(|name| !outputs.contains_key(*name)) {
            bail!("Study manifest is missing shared metadata or saved study objects");
        }
        for (relative, expected) in outputs {
            if !safe_relative(relative)
                || !(relative.starts_with("tables/") || relative.starts_with("plots/") || relative.starts_with("objects/"))
                || !expected.as_str().is_some_and(hash) {
                bail!("Invalid study output declaration: {relative}");
            }
            let mut output = directory.clone();
            for part in Path::new(relative).components() {
                output.push(part);
                if fs::symlink_metadata(&output)?.file_type().is_symlink() {
                    bail!("Symlink in study output: {}", output.display());
                }
            }
            regular_file(&output)?;
            if checksum(&output)? != expected.as_str().unwrap() {
                bail!("Study output checksum mismatch: {}", output.display());
            }
        }
        let effective = manifest["effective_dataset"].as_object().context("Study manifest missing effective_dataset")?;
        let original_batch = Path::new(effective.get("batch_dir").and_then(Value::as_str)
            .context("Study manifest has no original batch directory")?);
        if original_batch.file_name() != batch.file_name()
            || original_batch.parent().and_then(Path::file_name) != batch.parent().and_then(Path::file_name) {
            bail!("Study manifest belongs to another project or batch");
        }
        let (metadata_path, metadata) = checked_source(batch, original_batch,
            effective.get("metadata").and_then(Value::as_str).context("Study metadata path missing")?)?;
        let (quant_path, quant) = checked_source(batch, original_batch,
            effective.get("quantification").and_then(Value::as_str).context("Study quantification path missing")?)?;
        let inputs = manifest["inputs"].as_object().context("Study manifest missing input checksums")?;
        for (name, path) in [("metadata", &metadata_path), ("quantification", &quant_path)] {
            let expected = inputs.get(name).and_then(Value::as_str).filter(|value| hash(value)).context("Invalid study input SHA-256")?;
            if checksum(path)? != expected { bail!("Study input changed: {}", path.display()); }
        }
        let mut dataset = Value::Object(effective.clone());
        dataset["batch_dir"] = json!("../..");
        dataset["metadata"] = json!(metadata);
        dataset["quantification"] = json!(quant);
        if let Some(annotations) = effective.get("annotations") {
            let items = annotations.as_object().context("Invalid study annotations")?;
            let expected = inputs.get("annotations").context("Missing study annotation checksums")?;
            let mut relative = serde_json::Map::new();
            for (name, source) in items {
                if source.is_null() {
                    relative.insert(name.clone(), Value::Null);
                    continue;
                }
                let (path, name_relative) = checked_source(batch, original_batch,
                    source.as_str().context("Invalid study annotation path")?)?;
                let digest = expected.get(name).and_then(Value::as_str).filter(|value| hash(value)).context("Missing study annotation checksum")?;
                if checksum(&path)? != digest { bail!("Study annotation changed: {}", path.display()); }
                relative.insert(name.clone(), json!(name_relative));
            }
            dataset["annotations"] = Value::Object(relative);
        }
        dataset["input_sha256"] = json!({"metadata": inputs.get("metadata"), "quantification": inputs.get("quantification")});
        Ok(Self { hash: study_hash.to_owned(), dataset, metadata_path, quant_path })
    }
}
