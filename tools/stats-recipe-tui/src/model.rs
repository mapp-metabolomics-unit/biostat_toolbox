use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::fs::File;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

use crate::data::{self, Table};
use crate::study::Study;

pub struct Group {
    pub source: String,
    pub first: String,
    pub count: usize,
    pub included: bool,
    pub target: String,
}

pub struct Plan {
    pub folder: PathBuf,
    pub hash: String,
    pub existing: bool,
    pub dataset: Value,
    pub recipe: Value,
    pub corrected: Option<Table>,
    pub commands: Vec<String>,
}

impl Plan {
    fn verify_existing(&self) -> Result<()> {
        if !self.folder.is_dir() {
            bail!(
                "Hash-named recipe path is not a directory: {}",
                self.folder.display()
            );
        }
        let dataset_path = self.folder.join("dataset.yaml");
        let recipe_path = self.folder.join("recipe.yaml");
        let metadata_path = self.folder.join("metadata.tsv");
        for path in [&self.folder, &dataset_path, &recipe_path] {
            if fs::symlink_metadata(path)?.file_type().is_symlink() {
                bail!("Refusing symlink in existing recipe: {}", path.display());
            }
        }
        if self.corrected.is_some()
            && fs::symlink_metadata(&metadata_path)?
                .file_type()
                .is_symlink()
        {
            bail!(
                "Refusing symlink in existing recipe: {}",
                metadata_path.display()
            );
        }
        let dataset: Value = serde_yaml::from_reader(File::open(&dataset_path)?)?;
        let recipe: Value = serde_yaml::from_reader(File::open(&recipe_path)?)?;
        let metadata_matches = match &self.corrected {
            Some(table) => &Table::read_tsv(&metadata_path)? == table,
            None => !metadata_path.exists(),
        };
        if dataset != self.dataset || recipe != self.recipe || !metadata_matches {
            bail!("Hash-named recipe {} already exists with different contents; inspect it before choosing a different slug", self.folder.display());
        }
        Ok(())
    }

    pub fn save(&self) -> Result<()> {
        if self.existing {
            return self.verify_existing();
        }
        data::save(
            &self.folder,
            &self.dataset,
            &self.recipe,
            self.corrected.as_ref(),
        )
    }
}

pub struct Model {
    pub batch: PathBuf,
    pub study: Option<Study>,
    pub metadata_path: PathBuf,
    pub quant_path: PathBuf,
    pub metadata_paths: Vec<PathBuf>,
    pub quant_paths: Vec<PathBuf>,
    pub table: Table,
    pub aliases: [usize; 3],
    pub roles: [String; 3],
    pub measure: String,
    pub replacements: BTreeMap<String, String>,
    pub primary: usize,
    pub secondary: Option<usize>,
    pub display_group: Option<usize>,
    pub unit: Option<usize>,
    pub groups: Vec<Group>,
    pub slug: String,
    pub description: String,
    pub seed: u32,
    pub prep: Value,
    pub analyses: Value,
    pub annotation_candidates: BTreeMap<String, Vec<PathBuf>>,
    pub annotations: BTreeMap<String, Option<PathBuf>>,
    pub descriptive_only: bool,
    pub independent: bool,
    pub unpaired: bool,
    pub paired: bool,
}

fn valid_batch_name(name: &str, prefix: &str) -> bool {
    name.strip_prefix(prefix).is_some_and(|rest| {
        !rest.is_empty()
            && rest
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
    })
}

fn files(dir: &Path, extension: &str) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    if dir.is_dir() {
        for entry in fs::read_dir(dir).with_context(|| format!("Reading {}", dir.display()))? {
            let path = entry?.path();
            if path.is_file() && path.extension().is_some_and(|ext| ext == extension) {
                paths.push(path);
            }
        }
    }
    paths.sort();
    Ok(paths)
}

fn annotation_paths(batch: &Path) -> Result<BTreeMap<String, Vec<PathBuf>>> {
    let mut candidates = BTreeMap::new();
    for (key, path) in [
        (
            "canopus",
            batch.join("results/sirius/canopus_structure_summary.tsv"),
        ),
        (
            "structures",
            batch.join("results/sirius/structure_identifications.tsv"),
        ),
    ] {
        candidates.insert(
            key.to_owned(),
            if path.is_file() { vec![path] } else { vec![] },
        );
    }
    candidates.insert(
        "horizontal".to_owned(),
        files(&batch.join("results/tmp"), "tsv")?
            .into_iter()
            .filter(|path| {
                path.file_name().is_some_and(|n| {
                    n.to_string_lossy()
                        .ends_with("_met_annot_unified_horizontal.tsv")
                })
            })
            .collect(),
    );
    let mut enhancer = Vec::new();
    let root = batch.join("results/met_annot_enhancer");
    if root.is_dir() {
        for entry in fs::read_dir(&root).with_context(|| format!("Reading {}", root.display()))? {
            let path = entry?.path();
            if path.is_dir()
                && path
                    .file_name()
                    .is_some_and(|name| valid_batch_name(&name.to_string_lossy(), "mapp_batch_"))
            {
                enhancer.extend(files(&path, "tsv")?.into_iter().filter(|p| {
                    p.file_name().is_some_and(|n| {
                        n.to_string_lossy()
                            .ends_with("_spectral_match_results_repond_flat.tsv")
                    })
                }));
            }
        }
    }
    enhancer.sort();
    candidates.insert("met_annot_enhancer".to_owned(), enhancer);
    Ok(candidates)
}

fn words(name: &str) -> impl Iterator<Item = &str> {
    name.split(|c: char| !c.is_ascii_alphabetic())
        .filter(|word| !word.is_empty())
}

fn unit_name(name: &str) -> bool {
    words(name).any(|word| {
        [
            "id",
            "identifier",
            "subject",
            "specimen",
            "individual",
            "donor",
            "unit",
        ]
        .iter()
        .any(|w| word.eq_ignore_ascii_case(w))
    })
}

fn suggested_target(name: &str) -> String {
    let mut label = String::new();
    for ch in name.chars() {
        if ch.is_ascii_alphanumeric() || ch == '_' || ch == '-' {
            label.push(ch);
        } else if !label.ends_with('_') {
            label.push('_');
        }
    }
    let label = label.trim_matches('_');
    if valid_target(label) {
        label.to_owned()
    } else {
        "group".to_owned()
    }
}

fn valid_target(value: &str) -> bool {
    value != "NA"
        && value
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_alphabetic)
        && value
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-')
}

fn shell_quote(path: &Path) -> String {
    format!("'{}'", path.to_string_lossy().replace('\'', "'\\''"))
}

fn relative<'a>(batch: &Path, path: &'a Path) -> Result<&'a str> {
    path.strip_prefix(batch)
        .with_context(|| format!("{} is outside batch", path.display()))?
        .to_str()
        .context("Recipe paths must be valid UTF-8")
}

fn check_nonnegative(value: &Value, label: &str, max: Option<f64>) -> Result<()> {
    let number = value
        .as_f64()
        .with_context(|| format!("{label} must be a number"))?;
    if !number.is_finite() || number < 0.0 || max.is_some_and(|limit| number > limit) {
        bail!(
            "{label} must be a finite nonnegative number{}",
            max.map_or(String::new(), |n| format!(" <= {n}"))
        );
    }
    Ok(())
}

fn check_integer(value: &Value, label: &str, min: u64) -> Result<u64> {
    let number = value
        .as_u64()
        .with_context(|| format!("{label} must be a whole number"))?;
    if number < min || number > i32::MAX as u64 {
        bail!("{label} must be between {min} and {}", i32::MAX);
    }
    Ok(number)
}

fn check_method(section: &Value, name: &str, allowed: &[&str]) -> Result<()> {
    if !allowed
        .iter()
        .any(|method| section[name]["method"] == *method)
    {
        bail!("Invalid {name} method");
    }
    Ok(())
}

fn enabled(section: &Value, name: &str) -> Result<bool> {
    section[name]["enabled"]
        .as_bool()
        .with_context(|| format!("{name}.enabled must be a boolean"))
}

fn include_role(analyses: &Value, analysis: &str, option: &str) -> Result<bool> {
    match analyses[analysis].get(option) {
        None => Ok(false),
        Some(Value::Bool(value)) => Ok(*value),
        _ => bail!("{analysis}.{option} must be a boolean"),
    }
}

impl Model {
    pub fn load(batch: &Path) -> Result<Self> {
        Self::load_with_sources(batch, None, None)
    }

    fn load_with_sources(batch: &Path, metadata_source: Option<&Path>, quant_source: Option<&Path>) -> Result<Self> {
        let batch = batch
            .canonicalize()
            .with_context(|| format!("Batch {} not found", batch.display()))?;
        let parent = batch.parent().context("Batch has no project directory")?;
        if !batch.is_dir()
            || !batch
                .file_name()
                .and_then(|s| s.to_str())
                .is_some_and(|s| valid_batch_name(s, "mapp_batch_"))
            || !parent
                .file_name()
                .and_then(|s| s.to_str())
                .is_some_and(|s| valid_batch_name(s, "mapp_project_"))
        {
            bail!("Choose an existing mapp_project_*/mapp_batch_* directory.");
        }
        let metadata_paths = files(&batch.join("metadata/treated"), "tsv")?;
        let quant_paths = files(&batch.join("results/mzmine"), "csv")?;
        let metadata_path = metadata_source.map(Path::to_path_buf)
            .or_else(|| metadata_paths.first().cloned())
            .context("No treated metadata TSV files found")?;
        let quant_path = quant_source.map(Path::to_path_buf)
            .or_else(|| quant_paths.first().cloned())
            .context("No MZmine quantification CSV files found")?;
        let table = Table::read_tsv(&metadata_path)?;
        let annotation_candidates = annotation_paths(&batch)?;
        let annotations = annotation_candidates
            .iter()
            .map(|(key, paths)| (key.clone(), paths.first().cloned()))
            .collect();
        let mut model = Self {
            batch,
            metadata_path,
            quant_path,
            metadata_paths,
            quant_paths,
            table,
            aliases: [0; 3],
            roles: [String::new(), String::new(), String::new()],
            measure: String::new(),
            replacements: BTreeMap::new(),
            primary: 0,
            secondary: None,
            display_group: None,
            unit: None,
            groups: Vec::new(),
            study: None,
            slug: "batch-study".into(),
            description: String::new(),
            seed: 1,
            prep: json!({"blank_filter": {"enabled": true, "minimum_sample_to_blank_ratio": 5.0, "minimum_blank_prevalence": 0.5}, "qc_rsd_filter": {"enabled": true, "threshold_percent": 30.0}, "missing_values": {"method": "half_minimum"}, "normalization": {"method": "pqn"}, "transformation": {"method": "log2"}, "scaling": {"method": "pareto"}}),
            analyses: json!({"pca": {"enabled": true, "components": 5, "include_qc": false, "include_blanks": false}, "pcoa": {"enabled": true, "components": 5, "stage": "normalized", "distance": "bray", "include_qc": false, "include_blanks": false}, "omnibus": {"enabled": true}, "differential": {"enabled": true}, "plsda": {"enabled": false, "components": 3, "folds": 3, "repeats": 5, "permutations": 99}}),
            annotation_candidates,
            annotations,
            descriptive_only: false,
            independent: false,
            unpaired: false,
            paired: false,
        };
        model.reset_metadata_choices();
        model.reset_quant_choices()?;
        Ok(model)
    }

    pub fn load_comparison(batch: &Path, study_hash: &str) -> Result<Self> {
        let canonical = batch.canonicalize().with_context(|| format!("Batch {} not found", batch.display()))?;
        let study = Study::load(&canonical, study_hash)?;
        let mut model = Self::load_with_sources(&canonical, Some(&study.metadata_path), Some(&study.quant_path))?;
        model.reset_metadata_choices();
        for (slot, name) in ["filename", "sample_id", "sample_type"].iter().enumerate() {
            let column = study.dataset["metadata_columns"][name].as_str().unwrap_or(name);
            model.aliases[slot] = model.table.column(column)
                .with_context(|| format!("Study metadata lacks {column}"))?;
        }
        for (slot, name) in ["sample", "qc", "blank"].iter().enumerate() {
            model.roles[slot] = study.dataset["roles"]["levels"][name].as_str()
                .unwrap_or(["sample", "QC", "blank"][slot]).to_owned();
        }
        let fields = model.group_fields();
        if !fields.contains(&model.primary) {
            model.primary = fields.iter().copied()
                .find(|index| Some(*index) != model.unit)
                .or_else(|| fields.first().copied())
                .unwrap_or(0);
        }
        if model.unit.is_some_and(|index| !model.unit_fields().contains(&index)) {
            model.unit = None;
        }
        model.measure = study.dataset["intensity_measure"].as_str()
            .context("Study has no intensity measure")?.to_owned();
        model.replacements.clear();
        model.secondary = None;
        model.refresh_groups();
        model.slug = "biological-comparison".into();
        model.study = Some(study);
        Ok(model)
    }

    pub fn is_comparison(&self) -> bool { self.study.is_some() }

    fn reset_metadata_choices(&mut self) {
        // Reserve exact canonical columns before filling missing aliases from other fields.
        let canonical = ["filename", "sample_id", "sample_type"];
        let mut used: HashSet<usize> = canonical
            .iter()
            .filter_map(|name| self.table.column(name))
            .collect();
        for (slot, name) in canonical.iter().enumerate() {
            let index = self.table.column(name).unwrap_or_else(|| {
                (0..self.table.columns.len())
                    .find(|index| !used.contains(index))
                    .unwrap_or(0)
            });
            self.aliases[slot] = index;
            used.insert(index);
        }
        let levels = self.role_levels();
        let mut assigned = HashSet::new();
        for (index, preferred) in ["sample", "QC", "blank"].iter().enumerate() {
            let level = levels
                .iter()
                .find(|v| v.as_str() == *preferred && !assigned.contains(*v))
                .or_else(|| levels.iter().find(|v| !assigned.contains(*v)));
            self.roles[index] = level.cloned().unwrap_or_default();
            assigned.insert(self.roles[index].clone());
        }
        self.replacements.clear();
        self.secondary = None;
        self.display_group = None;
        let named_units: Vec<_> = self
            .unit_fields()
            .into_iter()
            .filter(|index| unit_name(&self.table.columns[*index]))
            .collect();
        self.unit = named_units.first().copied();
        let fields = self.group_fields();
        self.primary = fields
            .iter()
            .copied()
            .find(|i| {
                !named_units.contains(i)
                    && !words(&self.table.columns[*i]).any(|word| {
                        ["replicate", "technical", "batch"]
                            .iter()
                            .any(|w| word.eq_ignore_ascii_case(w))
                    })
            })
            .or_else(|| fields.iter().copied().find(|i| !named_units.contains(i)))
            .or_else(|| fields.first().copied())
            .unwrap_or(0);
        self.groups.clear();
        self.refresh_groups();
        self.independent = false;
        self.descriptive_only = false;
        self.paired = false;
        self.unpaired = false;
    }

    fn reset_quant_choices(&mut self) -> Result<()> {
        let header = data::quant_header(&self.quant_path)?;
        self.measure = ["height", "area"]
            .iter()
            .find(|measure| {
                header
                    .iter()
                    .any(|name| name.ends_with(&format!(" Peak {measure}")))
            })
            .context("Quantification has no Peak height or Peak area columns")?
            .to_string();
        self.replacements.clear();
        Ok(())
    }

    pub fn set_metadata(&mut self, index: usize) -> Result<()> {
        let path = self
            .metadata_paths
            .get(index)
            .context("Invalid metadata selection")?;
        let table = Table::read_tsv(path)?;
        self.metadata_path = path.clone();
        self.table = table;
        self.reset_metadata_choices();
        Ok(())
    }

    pub fn set_quant(&mut self, index: usize) -> Result<()> {
        let path = self
            .quant_paths
            .get(index)
            .context("Invalid quantification selection")?;
        let header = data::quant_header(path)?;
        let measure = ["height", "area"]
            .iter()
            .find(|measure| {
                header
                    .iter()
                    .any(|name| name.ends_with(&format!(" Peak {measure}")))
            })
            .context("Quantification has no Peak height or Peak area columns")?
            .to_string();
        self.quant_path = path.clone();
        self.measure = measure;
        self.replacements.clear();
        Ok(())
    }

    fn quant_names(&self) -> Result<BTreeSet<String>> {
        if !["height", "area"].contains(&self.measure.as_str()) {
            bail!("Choose Peak height or Peak area");
        }
        let suffix = format!(" Peak {}", self.measure);
        let mut names = BTreeSet::new();
        for name in data::quant_header(&self.quant_path)? {
            if let Some(sample) = name.strip_suffix(&suffix) {
                if sample.is_empty() || !names.insert(sample.to_owned()) {
                    bail!(
                        "Quantification sample names must be unique and nonempty for {}",
                        self.measure
                    );
                }
            }
        }
        if names.is_empty() {
            bail!("No Peak {} sample columns found", self.measure);
        }
        Ok(names)
    }

    pub fn alignment(&self) -> (Vec<String>, Vec<String>) {
        let metadata: BTreeSet<&str> = self
            .table
            .rows
            .iter()
            .filter_map(|row| row.get(self.aliases[0]).map(String::as_str))
            .collect();
        let quant = self.quant_names().unwrap_or_default();
        (
            metadata
                .iter()
                .filter(|name| !quant.contains(**name) && !self.replacements.contains_key(**name))
                .map(|s| (*s).to_owned())
                .collect(),
            quant
                .iter()
                .filter(|name| {
                    !metadata.contains(name.as_str())
                        && !self.replacements.values().any(|v| v == *name)
                })
                .cloned()
                .collect(),
        )
    }

    fn biological_rows(&self) -> impl Iterator<Item = &Vec<String>> {
        self.table.rows.iter().filter(|row| {
            row.get(self.aliases[2])
                .is_some_and(|role| role == &self.roles[0])
        })
    }

    fn complete_field(&self, index: usize) -> bool {
        index < self.table.columns.len()
            && self
                .biological_rows()
                .all(|row| !row[index].is_empty() && row[index] != "NA")
    }

    pub fn group_fields(&self) -> Vec<usize> {
        (0..self.table.columns.len())
            .filter(|index| {
                !self.aliases.contains(index)
                    && self.table.columns[*index] != "analysis_group"
                    && self.complete_field(*index)
                    && (2..=30).contains(
                        &self
                            .biological_rows()
                            .map(|row| &row[*index])
                            .collect::<HashSet<_>>()
                            .len(),
                    )
            })
            .collect()
    }

    pub fn second_fields(&self) -> Vec<usize> {
        (0..self.table.columns.len())
            .filter(|index| {
                *index != self.primary
                    && !self.aliases.contains(index)
                    && Some(*index) != self.unit
                    && self.table.columns[*index] != "analysis_group"
                    && self.complete_field(*index)
                    && self
                        .biological_rows()
                        .map(|row| &row[*index])
                        .collect::<HashSet<_>>()
                        .len()
                        <= 30
            })
            .collect()
    }

    pub fn unit_fields(&self) -> Vec<usize> {
        (0..self.table.columns.len())
            .filter(|index| {
                !self.aliases.contains(index)
                    && self.complete_field(*index)
                    && self.table.columns[*index] != "analysis_group"
            })
            .collect()
    }

    pub fn role_levels(&self) -> Vec<String> {
        self.table
            .rows
            .iter()
            .filter_map(|row| row.get(self.aliases[2]).cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    pub fn source_for_row(&self, row: &[String]) -> String {
        let first = row.get(self.primary).map_or("", String::as_str);
        match self.secondary.and_then(|index| row.get(index)) {
            Some(second) => {
                serde_json::to_string(&[first, second.as_str()]).expect("group names are strings")
            }
            None => first.to_owned(),
        }
    }

    pub fn refresh_groups(&mut self) {
        let old = std::mem::take(&mut self.groups)
            .into_iter()
            .map(|group| (group.source.clone(), group))
            .collect::<BTreeMap<_, _>>();
        let mut counts = BTreeMap::<String, (String, usize)>::new();
        if self.primary >= self.table.columns.len()
            || self
                .secondary
                .is_some_and(|i| i >= self.table.columns.len())
        {
            return;
        }
        for row in self.biological_rows() {
            let source = self.source_for_row(row);
            let group = counts
                .entry(source)
                .or_insert_with(|| (row[self.primary].clone(), 0));
            group.1 += 1;
        }
        let qualifying: HashSet<&str> = counts
            .values()
            .filter(|(_, count)| *count >= 3)
            .map(|(first, _)| first.as_str())
            .collect();
        let suggest_large = qualifying.len() >= 2 && counts.values().any(|(_, count)| *count < 3);
        self.groups = counts
            .into_iter()
            .map(|(source, (first, count))| {
                if let Some(previous) = old.get(&source) {
                    Group {
                        source,
                        first,
                        count,
                        included: previous.included,
                        target: previous.target.clone(),
                    }
                } else {
                    Group {
                        source,
                        target: suggested_target(&first),
                        first,
                        count,
                        included: !suggest_large || count >= 3,
                    }
                }
            })
            .collect();
        self.independent = false;
        self.paired = false;
        self.unpaired = false;
    }

    pub fn dependency_warnings(&self) -> Vec<String> {
        let Some(unit) = self.unit else {
            return vec![
                "No individual-unit ID column chosen; paired injections cannot be detected.".into(),
            ];
        };
        if unit >= self.table.columns.len() {
            return vec!["Select a valid experimental-unit ID column.".into()];
        }
        let included: HashSet<&str> = self
            .groups
            .iter()
            .filter(|g| g.included)
            .map(|g| g.source.as_str())
            .collect();
        let mut ids = BTreeMap::<&str, usize>::new();
        let mut prefixes = BTreeMap::<&str, BTreeSet<&str>>::new();
        for row in self
            .biological_rows()
            .filter(|row| included.contains(self.source_for_row(row).as_str()))
        {
            let id = row[unit].as_str();
            *ids.entry(id).or_default() += 1;
            if let Some(split) = id.rfind(|c| c == '-' || c == '_') {
                let (prefix, suffix) = id.split_at(split);
                let suffix = &suffix[1..];
                if prefix.as_bytes().last().is_some_and(u8::is_ascii_digit)
                    && suffix
                        .as_bytes()
                        .first()
                        .is_some_and(u8::is_ascii_alphabetic)
                    && suffix.bytes().all(|b| b.is_ascii_alphanumeric())
                {
                    prefixes.entry(prefix).or_default().insert(id);
                }
            }
        }
        let repeated: Vec<_> = ids
            .iter()
            .filter(|(_, count)| **count > 1)
            .map(|(id, _)| *id)
            .collect();
        let mut warnings = Vec::new();
        if !repeated.is_empty() {
            warnings.push(format!(
                "POSSIBLE DEPENDENCY: repeated {} values: {}",
                self.table.columns[unit],
                repeated.join(", ")
            ));
        }
        for (prefix, distinct) in prefixes
            .into_iter()
            .filter(|(_, distinct)| distinct.len() > 1)
        {
            warnings.push(format!("POSSIBLE DEPENDENCY: {} prefix {prefix} appears in {}; verify whether these are one specimen.", self.table.columns[unit], distinct.into_iter().collect::<Vec<_>>().join(", ")));
        }
        warnings
    }

    pub fn assign_filename(&mut self, old: &str, new: &str, verification: &str) -> Result<()> {
        if verification != "VERIFY" {
            bail!("Type VERIFY to confirm the exact injection correspondence");
        }
        let (metadata_only, quant_only) = self.alignment();
        if !metadata_only.iter().any(|v| v == old) || !quant_only.iter().any(|v| v == new) {
            bail!("Choose a currently unmatched metadata filename and an unused quantification filename");
        }
        self.replacements.insert(old.to_owned(), new.to_owned());
        Ok(())
    }

    fn verify_unblocked_units(&self, allow_shared: bool) -> Result<()> {
        let unit = self.unit.context("Choose a verified biological-unit ID column to rule out repeated injections within groups")?;
        if !self.unit_fields().contains(&unit) {
            bail!("Choose a complete biological-unit ID column");
        }
        if unit == self.primary || self.secondary == Some(unit) {
            bail!("Biological-unit ID column must not be a grouping column");
        }
        let selected = self.groups.iter().filter(|group| group.included)
            .map(|group| (group.source.as_str(), group.target.as_str()))
            .collect::<BTreeMap<_, _>>();
        if allow_shared && selected.values().copied().collect::<BTreeSet<_>>().len() != 2 {
            bail!("UNPAIRED mode requires exactly two selected target groups");
        }
        let mut seen = HashSet::<(&str, &str)>::new();
        let mut first_group = BTreeMap::<&str, &str>::new();
        for row in self.biological_rows() {
            let source = self.source_for_row(row);
            let Some(&target) = selected.get(source.as_str()) else { continue };
            let id = row[unit].as_str();
            if id.is_empty() || id == "NA" {
                bail!("Biological-unit ID is missing in selected group {target}");
            }
            if !seen.insert((id, target)) {
                bail!("Biological unit {id} has repeated injections within group {target}; use one biological observation per group");
            }
            if let Some(previous) = first_group.insert(id, target) {
                if !allow_shared && previous != target {
                    bail!("Biological unit {id} occurs in groups {previous} and {target}; choose UNPAIRED to ignore matching or PAIRED to model it, not INDEPENDENT");
                }
            }
        }
        Ok(())
    }

    pub fn confirm_unpaired(&mut self, token: &str) -> Result<()> {
        if token != "UNPAIRED" {
            bail!("Type UNPAIRED to ignore matching only after verifying one biological observation per group");
        }
        self.verify_unblocked_units(true)?;
        self.unpaired = true;
        self.independent = false;
        self.paired = false;
        self.descriptive_only = false;
        Ok(())
    }

    fn verify_pairs(&self) -> Result<()> {
        let unit = self.unit.context("Paired comparison requires a verified biological-unit ID column")?;
        if !self.unit_fields().contains(&unit) {
            bail!("Choose a complete biological-unit ID column for pairing");
        }
        if unit == self.primary || self.secondary == Some(unit) {
            bail!("Pairing column must identify biological units, not a grouping column");
        }
        let targets = self.groups.iter().filter(|group| group.included)
            .map(|group| group.target.as_str()).collect::<BTreeSet<_>>();
        if targets.len() != 2 {
            bail!("Paired comparison requires exactly two selected target groups");
        }
        let labels: Vec<_> = targets.into_iter().collect();
        let group_targets = self.groups.iter().filter(|group| group.included)
            .map(|group| (group.source.as_str(), group.target.as_str()))
            .collect::<BTreeMap<_, _>>();
        let mut blocks = BTreeMap::<&str, [usize; 2]>::new();
        for row in self.biological_rows() {
            let source = self.source_for_row(row);
            if let Some(target) = group_targets.get(source.as_str()) {
                let block = row[unit].as_str();
                if block.is_empty() || block == "NA" {
                    bail!("Paired comparison has a missing biological-unit ID in selected group {target}");
                }
                let index = usize::from(*target == labels[1]);
                let counts = blocks.entry(block).or_default();
                counts[index] += 1;
                if counts[index] > 1 {
                    bail!("Paired biological unit {block} has multiple samples in group {target}; exactly one per group is required");
                }
            }
        }
        for (block, counts) in &blocks {
            if counts != &[1, 1] {
                bail!("Paired biological unit {block} is missing a sample in group {}; every selected unit needs one sample per group", labels[usize::from(counts[0] == 1)]);
            }
        }
        if blocks.len() < 3 {
            bail!("Paired comparison requires at least three complete matched biological units");
        }
        Ok(())
    }

    pub fn confirm_paired(&mut self, token: &str) -> Result<()> {
        if token != "PAIRED" {
            bail!("Type PAIRED only after verifying that the selected column identifies the same biological unit in both groups");
        }
        self.verify_pairs()?;
        self.paired = true;
        self.independent = false;
        self.unpaired = false;
        self.descriptive_only = false;
        self.analyses["plsda"]["enabled"] = json!(false);
        Ok(())
    }

    pub fn confirm_independence(&mut self, token: &str) -> Result<()> {
        if token != "INDEPENDENT" {
            bail!("Type INDEPENDENT only after verifying independent experimental units");
        }
        if self.unit.is_some() {
            self.verify_unblocked_units(false)?;
        }
        self.paired = false;
        self.unpaired = false;
        self.independent = true;
        self.descriptive_only = false;
        Ok(())
    }

    fn validate_settings(&self, target_counts: &BTreeMap<&str, usize>) -> Result<()> {
        let analyses = &self.analyses;
        if self.is_comparison() {
            for key in ["omnibus", "differential", "plsda"] {
                enabled(analyses, key)?;
            }
            if self.paired && enabled(analyses, "plsda")? {
                bail!("PLS-DA is not supported for paired comparisons; disable PLS-DA.");
            }
            check_integer(&analyses["plsda"]["components"], "PLS-DA components", 1)?;
            let folds = check_integer(&analyses["plsda"]["folds"], "PLS-DA folds", 2)?;
            check_integer(&analyses["plsda"]["repeats"], "PLS-DA repeats", 2)?;
            check_integer(&analyses["plsda"]["permutations"], "PLS-DA permutations", 19)?;
            if enabled(analyses, "plsda")?
                && (target_counts.values().copied().min().unwrap_or(0) < 4
                    || folds > target_counts.values().copied().min().unwrap_or(0) as u64)
            {
                bail!("PLS-DA needs >=4 replicates per group and folds <= the smallest group.");
            }
            if self.descriptive_only
                && ["omnibus", "differential", "plsda"].iter()
                    .any(|key| analyses[*key]["enabled"] == true)
            {
                bail!("Descriptive-only plans cannot enable inferential analyses.");
            }
            if !self.independent && !self.unpaired && !self.paired && !self.descriptive_only {
                bail!("Choose descriptive-only, confirm UNPAIRED or PAIRED, or type INDEPENDENT after verifying study design");
            }
            return Ok(());
        }
        let prep = &self.prep;
        let blank = enabled(prep, "blank_filter")?;
        let qc = enabled(prep, "qc_rsd_filter")?;
        check_nonnegative(&prep["blank_filter"]["minimum_sample_to_blank_ratio"], "Minimum sample/blank ratio", None)?;
        check_nonnegative(&prep["blank_filter"]["minimum_blank_prevalence"], "Blank prevalence", Some(1.0))?;
        check_nonnegative(&prep["qc_rsd_filter"]["threshold_percent"], "QC RSD threshold", None)?;
        for (key, values) in [
            ("missing_values", &["none", "half_minimum"][..]),
            ("normalization", &["none", "total_sum", "median", "pqn"]),
            ("transformation", &["none", "log2", "log10", "log1p"]),
            ("scaling", &["none", "pareto", "autoscale"]),
        ] {
            check_method(prep, key, values)?;
        }
        let qc_count = self.biological_role_count(1);
        let blank_count = self.biological_role_count(2);
        if blank && blank_count == 0 {
            bail!("Blank filter needs at least one blank; disable it in settings.");
        }
        if qc && qc_count < 3 {
            bail!("QC RSD filter needs at least three QC samples; disable it in settings.");
        }
        for analysis in ["pca", "pcoa"] {
            enabled(analyses, analysis)?;
            let include_qc = include_role(analyses, analysis, "include_qc")?;
            let include_blanks = include_role(analyses, analysis, "include_blanks")?;
            if enabled(analyses, analysis)? {
                if include_qc && qc_count == 0 {
                    bail!("{analysis}.include_qc needs at least one QC injection.");
                }
                if include_blanks && blank_count == 0 {
                    bail!("{analysis}.include_blanks needs at least one blank injection.");
                }
            }
        }
        check_integer(&analyses["pca"]["components"], "PCA components", 1)?;
        check_integer(&analyses["pcoa"]["components"], "PCoA components", 1)?;
        if !["imputed", "normalized", "transformed", "scaled"]
            .iter().any(|v| analyses["pcoa"]["stage"] == *v)
            || !["bray", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"]
                .iter().any(|v| analyses["pcoa"]["distance"] == *v)
        {
            bail!("Choose a supported PCoA stage and distance");
        }
        if enabled(analyses, "pcoa")?
            && analyses["pcoa"]["distance"] == "bray"
            && analyses["pcoa"]["stage"] != "normalized"
            && analyses["pcoa"]["stage"] != "imputed"
        {
            bail!("Bray-Curtis requires a nonnegative imputed/normalized PCoA stage.");
        }
        Ok(())
    }

    fn biological_role_count(&self, slot: usize) -> usize {
        self.table
            .rows
            .iter()
            .filter(|row| row[self.aliases[2]] == self.roles[slot])
            .count()
    }

    pub fn plan(&self) -> Result<Plan> {
        let columns = &self.table.columns;
        if self.aliases.iter().any(|i| *i >= columns.len())
            || self.aliases.iter().collect::<HashSet<_>>().len() != 3
            || ["filename", "sample_id", "sample_type"]
                .iter()
                .enumerate()
                .any(|(slot, canonical)| {
                    columns
                        .iter()
                        .position(|c| c == canonical)
                        .is_some_and(|i| i != self.aliases[slot])
                })
        {
            bail!("Metadata aliases must be distinct and cannot collide with canonical columns.");
        }
        for (slot, name) in [(0, "filename"), (1, "sample_id")] {
            let mut seen = HashSet::new();
            if self.table.rows.iter().any(|row| {
                row.get(self.aliases[slot])
                    .is_none_or(|v| v.is_empty() || v == "NA" || !seen.insert(v))
            }) {
                bail!("Metadata {name} must be unique and nonempty.");
            }
        }
        let levels = self.role_levels();
        if self.roles.iter().any(String::is_empty)
            || self.roles.iter().collect::<HashSet<_>>().len() != 3
            || levels.iter().any(|level| !self.roles.contains(level))
        {
            bail!(
                "Sample roles must be distinct and every metadata row must have a recognized role."
            );
        }
        let quant_names = self.quant_names()?;
        let metadata_names: BTreeSet<&str> = self
            .table
            .rows
            .iter()
            .map(|r| r[self.aliases[0]].as_str())
            .collect();
        for (old, new) in &self.replacements {
            if !metadata_names.contains(old.as_str())
                || quant_names.contains(old)
                || metadata_names.contains(new.as_str())
                || !quant_names.contains(new)
            {
                bail!("Filename correction {old} -> {new} is not an exact unmatched pair");
            }
        }
        let translated: BTreeSet<&str> = self
            .table
            .rows
            .iter()
            .map(|row| {
                self.replacements
                    .get(&row[self.aliases[0]])
                    .map_or(row[self.aliases[0]].as_str(), String::as_str)
            })
            .collect();
        if translated != quant_names.iter().map(String::as_str).collect()
            || translated.len() != self.table.rows.len()
        {
            bail!("Unmatched injections remain: align exact filenames and VERIFY each correction before saving");
        }
        if columns.iter().any(|name| name == "analysis_group") {
            bail!("Original metadata already contains analysis_group");
        }
        if self.is_comparison()
            && self.secondary.is_some_and(|index| !self.second_fields().contains(&index))
        {
            bail!("Choose a complete secondary group column");
        }
        if self.is_comparison() && !self.group_fields().contains(&self.primary) {
            bail!("Choose a complete categorical primary group column (2–30 levels)");
        }
        if self.is_comparison() && self
            .unit
            .is_some_and(|index| !self.unit_fields().contains(&index))
        {
            bail!("Choose a complete experimental-unit ID column");
        }
        if self.biological_role_count(0) == 0 {
            bail!("A study requires at least one biological sample");
        }
        let mut selected = Vec::new();
        let mut target_counts = BTreeMap::<&str, usize>::new();
        if self.is_comparison() {
            let mut actual = BTreeMap::<String, (String, usize)>::new();
            for row in self.biological_rows() {
                let group = actual.entry(self.source_for_row(row))
                    .or_insert_with(|| (row[self.primary].clone(), 0));
                group.1 += 1;
            }
            if self.groups.len() != actual.len() {
                bail!("Biological groups changed; refresh group choices");
            }
            let mut seen = HashSet::new();
            for group in &self.groups {
                if !seen.insert(group.source.as_str())
                    || actual.get(&group.source) != Some(&(group.first.clone(), group.count))
                {
                    bail!("Biological groups changed; refresh group choices");
                }
                if group.included {
                    if !valid_target(&group.target) {
                        bail!("Target label must start with a letter and use letters, digits, _ or -, not NA");
                    }
                    *target_counts.entry(&group.target).or_default() += group.count;
                    selected.push(json!({"source": group.source, "target": group.target}));
                }
            }
            if target_counts.len() < 2 || target_counts.values().any(|count| *count < 3) {
                bail!("At least two target groups with at least three biological samples each are required.");
            }
        }
        if self.paired {
            if !self.is_comparison() {
                bail!("Paired design is only available for comparisons");
            }
            self.verify_pairs()?;
        }
        if self.unpaired || (self.independent && self.unit.is_some()) {
            self.verify_unblocked_units(self.unpaired)?;
        }
        self.validate_settings(&target_counts)?;
        let slug = self.slug.as_bytes();
        if !(3..=51).contains(&slug.len())
            || !slug.first().is_some_and(u8::is_ascii_lowercase)
            || !slug
                .iter()
                .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || *b == b'_' || *b == b'-')
        {
            bail!("Recipe directory must be 3–51 lowercase letters/digits/_/-, starting with a letter");
        }
        if self.description.chars().count() > 160 || self.description.chars().any(char::is_control)
        {
            bail!("Description must be one line of at most 160 characters");
        }
        if self.seed > i32::MAX as u32 {
            bail!("Seed must be <= {}", i32::MAX);
        }
        // Hash authored settings and any verified metadata corrections before embedding the folder name.
        let copied = !self.is_comparison() && !self.replacements.is_empty();
        let mut dataset = if let Some(study) = &self.study {
            study.dataset.clone()
        } else {
            let parent = self.batch.parent().context("Batch has no project")?;
            let mut dataset = json!({"id": format!("{}/{}", parent.file_name().unwrap_or_default().to_string_lossy(), self.batch.file_name().unwrap_or_default().to_string_lossy()),
                "batch_dir": "../..", "metadata": if copied { "metadata.tsv".to_owned() } else { relative(&self.batch, &self.metadata_path)?.to_owned() },
                "quantification": relative(&self.batch, &self.quant_path)?, "feature_id_column": "row ID", "intensity_measure": self.measure});
            let nonstandard = ["filename", "sample_id", "sample_type"]
                .iter()
                .enumerate()
                .filter(|(index, name)| **name != columns[self.aliases[*index]])
                .map(|(index, name)| ((*name).to_owned(), json!(columns[self.aliases[index]])))
                .collect::<serde_json::Map<String, Value>>();
            if !nonstandard.is_empty() {
                dataset["metadata_columns"] = Value::Object(nonstandard);
            }
            if self.roles != ["sample", "QC", "blank"] {
                dataset["roles"] = json!({"column": "sample_type", "levels": {"sample": self.roles[0], "qc": self.roles[1], "blank": self.roles[2]}});
            }
            let mut annotations = serde_json::Map::new();
            for (key, path) in &self.annotations {
                let choices = self.annotation_candidates.get(key).context("Unknown annotation type")?;
                if let Some(path) = path {
                    if !choices.contains(path) || !path.is_file() {
                        bail!("Annotation file {} is unavailable", path.display());
                    }
                    annotations.insert(key.clone(), json!(relative(&self.batch, path)?));
                }
            }
            if !annotations.is_empty() {
                dataset["annotations"] = Value::Object(annotations);
            }
            dataset
        };
        let mut recipe = if let Some(study) = &self.study {
            let mut sample_plan = json!({"group_column": columns[self.primary], "groups": selected});
            if let Some(index) = self.secondary {
                sample_plan["second_column"] = json!(columns[index]);
            }
            let mut comparison = json!({"kind": "comparison", "study_hash": study.hash, "seed": self.seed,
                "sample_plan": sample_plan,
                "design": {"group": "analysis_group", "contrasts": "all"},
                "analyses": {"omnibus": self.analyses["omnibus"], "differential": self.analyses["differential"], "plsda": self.analyses["plsda"]}});
            if self.paired {
                comparison["design"]["block"] = json!(columns[self.unit.expect("verified by verify_pairs")]);
                comparison["design"]["block_verified"] = json!(true);
            }
            comparison
        } else {
            json!({"kind": "study", "seed": self.seed, "preprocessing": self.prep,
                "analyses": {"pca": self.analyses["pca"], "pcoa": self.analyses["pcoa"]}})
        };
        if !self.is_comparison() {
            if let Some(index) = self.display_group {
                if !self.group_fields().contains(&index) {
                    bail!("Choose a complete biological display-group column");
                }
                recipe["display_group"] = json!(columns[index]);
            }
        }
        if !self.description.is_empty() {
            recipe["presentation"] = json!({"description": self.description});
        }
        let corrected = if copied {
            let mut table = self.table.clone();
            for row in &mut table.rows {
                if let Some(replacement) = self.replacements.get(&row[self.aliases[0]]) {
                    row[self.aliases[0]] = replacement.clone();
                }
            }
            Some(table)
        } else {
            None
        };
        let fingerprint = serde_json::to_vec(&(
            &dataset,
            &recipe,
            corrected
                .as_ref()
                .map(|table| (&table.columns, &table.rows)),
        ))?;
        let hash = format!("{:x}", Sha256::digest(&fingerprint));
        let name = format!("{}-{}", self.slug, &hash[..12]);
        let folder = self.batch.join("recipes").join(&name);
        if copied {
            dataset["metadata"] = json!(format!("recipes/{name}/metadata.tsv"));
        }
        let biostats = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../biostats");
        let prefix = format!("cd {} && rig run -r 4.6.1 -f run.R -- ", shell_quote(&biostats));
        let options = format!(
            " --dataset {} --recipe {}",
            shell_quote(&folder.join("dataset.yaml")),
            shell_quote(&folder.join("recipe.yaml"))
        );
        let commands = ["validate", "plan", "run"]
            .iter()
            .map(|action| format!("{prefix}{action}{options}"))
            .collect::<Vec<_>>();
        let plan = Plan {
            folder,
            hash,
            existing: false,
            dataset,
            recipe,
            corrected,
            commands,
        };
        let existing = fs::symlink_metadata(&plan.folder)
            .map(|_| true)
            .or_else(|error| {
                if error.kind() == std::io::ErrorKind::NotFound {
                    Ok(false)
                } else {
                    Err(error)
                }
            })?;
        if existing {
            plan.verify_existing()?;
        }
        Ok(Plan { existing, ..plan })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn study_and_incremental_comparisons_reuse_corrected_shared_metadata() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let batch = temp.path().join("mapp_project_test/mapp_batch_test");
        fs::create_dir_all(batch.join("metadata/treated"))?;
        fs::create_dir_all(batch.join("results/mzmine"))?;
        let source = batch.join("metadata/treated/metadata.tsv");
        fs::write(&source, concat!(
            "filename\tsample_id\tsample_type\tindividual_id\tcondition\tphase\tcollection_date\n",
            "A1.mzML\tA1\tsample\tA1\tcase\tall\t2026-01-01\nA2.mzML\tA2\tsample\tA2\tcase\tall\t2026-01-02\nA3.mzML\tA3\tsample\tA3\tcase\tall\t2026-01-03\nA4.mzML\tA4\tsample\tA4\tcase\tall\t2026-01-04\n",
            "B1.mzML\tB1\tsample\tA1\tcontrol\tall\t2026-01-01\nB2.mzML\tB2\tsample\tA2\tcontrol\tall\t2026-01-02\nB3.mzML\tB3\tsample\tA3\tcontrol\tall\t2026-01-03\nB4.mzML\tB4\tsample\tA4\tcontrol\tall\t2026-01-04\n",
            "Q1.mzML\tQ1\tQC\tQ1\tNA\tNA\tNA\nQ2.mzML\tQ2\tQC\tQ2\tNA\tNA\tNA\nQ3.mzML\tQ3\tQC\tQ3\tNA\tNA\tNA\n",
            "blank_old.mzML\tblank\tblank\tblank\tNA\tNA\tNA\n"
        ))?;
        fs::write(batch.join("results/mzmine/quant.csv"), concat!(
            "row ID,A1.mzML Peak height,A2.mzML Peak height,A3.mzML Peak height,A4.mzML Peak height,",
            "B1.mzML Peak height,B2.mzML Peak height,B3.mzML Peak height,B4.mzML Peak height,",
            "Q1.mzML Peak height,Q2.mzML Peak height,Q3.mzML Peak height,blank_new.mzML Peak height,\n"
        ))?;
        let mut model = Model::load(&batch)?;
        model.assign_filename("blank_old.mzML", "blank_new.mzML", "VERIFY")?;
        let plan = model.plan()?; // Study does not require inferential-group certification.
        assert_eq!(plan.recipe["kind"], "study");
        assert!(plan.recipe.get("sample_plan").is_none());
        assert!(plan.recipe.get("design").is_none());
        assert_eq!(plan.recipe["analyses"].as_object().unwrap().len(), 2);
        assert!(plan.recipe.get("display_group").is_none());
        let mut one_group = Model::load(&batch)?;
        one_group.assign_filename("blank_old.mzML", "blank_new.mzML", "VERIFY")?;
        for row in &mut one_group.table.rows {
            if row[2] == "sample" { row[4] = "only_group".into(); }
        }
        assert_eq!(one_group.plan()?.recipe["kind"], "study");
        plan.save()?;
        let corrected_path = plan.folder.join("metadata.tsv");
        let corrected = Table::read_tsv(&corrected_path)?;
        assert_eq!(corrected.rows.last().unwrap()[0], "blank_new.mzML");
        assert_eq!(Table::read_tsv(&source)?.rows.last().unwrap()[0], "blank_old.mzML");
        assert!(model.plan()?.existing);

        let study_hash = "a".repeat(64);
        let study_dir = batch.join("results/stats").join(&study_hash);
        fs::create_dir_all(study_dir.join("tables"))?;
        fs::create_dir_all(study_dir.join("objects"))?;
        fs::write(study_dir.join("objects/dataset.rds"), "test dataset object")?;
        fs::write(study_dir.join("objects/processed.rds"), "test processed object")?;
        fs::copy(&corrected_path, study_dir.join("tables/sample_metadata.tsv"))?;
        let quant = batch.join("results/mzmine/quant.csv").canonicalize()?;
        let metadata = corrected_path.canonicalize()?;
        let digest = |path: &Path| -> Result<String> {
            Ok(format!("{:x}", Sha256::digest(fs::read(path)?)))
        };
        let output = study_dir.join("tables/sample_metadata.tsv");
        let mut effective_dataset = plan.dataset.clone();
        effective_dataset["batch_dir"] = json!(batch.canonicalize()?.to_str());
        effective_dataset["metadata"] = json!(metadata.to_str());
        effective_dataset["quantification"] = json!(quant.to_str());
        let manifest = json!({
            "kind": "study", "schema_version": "1.0.0", "status": "complete",
            "hashes": {"run_hash": study_hash, "dataset_hash": "b".repeat(64), "preprocess_hash": "c".repeat(64)},
            "inputs": {"metadata": digest(&metadata)?, "quantification": digest(&quant)?, "annotations": {}},
            "effective_dataset": effective_dataset,
            "effective_recipe": {"kind": "study", "preprocessing": plan.recipe["preprocessing"], "analyses": plan.recipe["analyses"]},
            "outputs": {"tables/sample_metadata.tsv": digest(&output)?,
                        "objects/dataset.rds": digest(&study_dir.join("objects/dataset.rds"))?,
                        "objects/processed.rds": digest(&study_dir.join("objects/processed.rds"))?}
        });
        fs::write(study_dir.join("manifest.yaml"), serde_json::to_vec(&manifest)?)?;
        fs::write(study_dir.join("COMPLETE"), serde_json::to_vec(&json!({
            "run_hash": study_hash, "content_sha256": "d".repeat(64)
        }))?)?;
        assert!(Model::load_comparison(&batch, "a").is_err());
        let mut comparison = Model::load_comparison(&batch, &study_hash)?;
        assert_eq!(comparison.metadata_path, metadata);
        assert!(comparison.plan().is_err()); // Inferential designs require explicit certification.
        let selected_unit = comparison.unit.take();
        assert!(comparison.confirm_unpaired("UNPAIRED").err().unwrap().to_string().contains("biological-unit ID"));
        comparison.confirm_independence("INDEPENDENT")?;
        assert!(comparison.plan()?.recipe["design"].get("block").is_none());
        comparison.unit = selected_unit;
        assert!(comparison.plan().err().unwrap().to_string().contains("UNPAIRED"));
        assert!(comparison.confirm_independence("INDEPENDENT").err().unwrap().to_string().contains("UNPAIRED"));
        assert!(comparison.confirm_unpaired("INDEPENDENT").is_err());
        comparison.table.rows[1][3] = "A1".into();
        assert!(comparison.confirm_unpaired("UNPAIRED").err().unwrap().to_string().contains("repeated injections within group"));
        comparison.table.rows[1][3] = "A2".into();
        for index in 4..8 {
            comparison.table.rows[index][3] = format!("B{}", index - 3);
        }
        comparison.confirm_independence("INDEPENDENT")?;
        for index in 4..8 {
            comparison.table.rows[index][3] = format!("A{}", index - 3);
        }
        comparison.confirm_unpaired("UNPAIRED")?;
        let first = comparison.plan()?;
        assert_eq!(first.recipe["kind"], "comparison");
        assert_eq!(first.recipe["study_hash"], study_hash);
        assert_eq!(first.recipe["sample_plan"]["groups"].as_array().unwrap().len(), 2);
        assert!(first.recipe.get("preprocessing").is_none());
        assert!(first.recipe["analyses"].get("pca").is_none());
        assert!(first.recipe["analyses"].get("pcoa").is_none());
        assert!(first.corrected.is_none());
        assert_eq!(first.dataset["metadata"], plan.dataset["metadata"]);
        assert_eq!(first.dataset["input_sha256"]["metadata"], digest(&metadata)?);
        first.save()?;
        assert!(first.recipe["design"].get("block").is_none());
        assert!(comparison.unpaired);
        assert!(!comparison.independent);
        assert!(comparison.confirm_paired("INDEPENDENT").is_err());
        comparison.confirm_paired("PAIRED")?;
        let matched = comparison.plan()?;
        assert_eq!(matched.recipe["design"]["block"], "individual_id");
        assert_eq!(matched.recipe["design"]["block_verified"], true);
        assert_ne!(first.hash, matched.hash);
        assert!(!comparison.independent);
        assert_eq!(matched.recipe["analyses"]["plsda"]["enabled"], false);
        comparison.unit = comparison.table.column("collection_date");
        comparison.confirm_paired("PAIRED")?;
        assert_eq!(comparison.plan()?.recipe["design"]["block"], "collection_date");
        comparison.table.rows[4][6] = "2026-01-02".into();
        assert!(comparison.plan().err().unwrap().to_string().contains("multiple samples"));
        comparison.table.rows[4][6] = "2026-01-01".into();
        comparison.unit = comparison.table.column("individual_id");
        comparison.analyses["plsda"]["enabled"] = json!(true);
        assert!(comparison.plan().err().unwrap().to_string().contains("PLS-DA"));
        comparison.analyses["plsda"]["enabled"] = json!(false);
        comparison.table.rows[4][3] = "A2".into();
        assert!(comparison.plan().err().unwrap().to_string().contains("multiple samples"));
        assert!(comparison.confirm_paired("PAIRED").err().unwrap().to_string().contains("multiple samples"));
        comparison.table.rows[4][3] = "unmatched".into();
        assert!(comparison.plan().err().unwrap().to_string().contains("missing a sample"));
        comparison.table.rows[4][3].clear();
        assert!(comparison.plan().err().unwrap().to_string().contains("complete experimental-unit ID"));
        comparison.table.rows[4][3] = "A1".into();
        comparison.unit = None;
        assert!(comparison.plan().err().unwrap().to_string().contains("biological-unit ID"));
        comparison.unit = comparison.table.column("individual_id");
        comparison.primary = comparison.unit.unwrap();
        assert!(comparison.verify_pairs().err().unwrap().to_string().contains("Pairing column must identify"));
        comparison.primary = comparison.table.column("condition").unwrap();
        let rows = comparison.table.rows.clone();
        comparison.table.rows.retain(|row| !["A4", "B4"].contains(&row[1].as_str()));
        comparison.refresh_groups();
        comparison.confirm_paired("PAIRED")?; // Three complete units are sufficient.
        comparison.table.rows.retain(|row| !["A3", "B3"].contains(&row[1].as_str()));
        comparison.refresh_groups();
        assert!(comparison.confirm_paired("PAIRED").unwrap_err().to_string().contains("at least three complete"));
        comparison.table.rows = rows;
        comparison.refresh_groups();
        comparison.confirm_unpaired("UNPAIRED")?;
        assert!(comparison.plan()?.recipe["design"].get("block").is_none());
        comparison.secondary = comparison.table.column("phase");
        comparison.refresh_groups();
        comparison.confirm_unpaired("UNPAIRED")?;
        let secondary = comparison.plan()?;
        assert_eq!(secondary.recipe["sample_plan"]["second_column"], "phase");
        assert_eq!(secondary.recipe["sample_plan"]["groups"][0]["source"], "[\"case\",\"all\"]");
        assert!(secondary.corrected.is_none());
        assert_eq!(secondary.dataset["metadata"], first.dataset["metadata"]);
        assert_ne!(first.hash, secondary.hash);
        comparison.groups[0].target = "renamed".into();
        let second = comparison.plan()?;
        assert_ne!(first.hash, second.hash);
        assert_eq!(fs::read(&output)?, fs::read(&corrected_path)?);
        let relocated = temp.path().join("new-checkout/mapp_project_test/mapp_batch_test");
        fs::create_dir_all(relocated.parent().unwrap())?;
        fs::rename(&batch, &relocated)?;
        let moved = Model::load_comparison(&relocated, &study_hash)?;
        assert_eq!(moved.metadata_path, relocated.join("recipes")
            .join(plan.folder.file_name().unwrap()).join("metadata.tsv"));
        fs::write(relocated.join("results/stats").join(&study_hash)
            .join("tables/sample_metadata.tsv"), "tampered")?;
        assert!(Model::load_comparison(&relocated, &study_hash).is_err());
        Ok(())
    }
}
