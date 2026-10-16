use std::io::{self, Stdout};
use std::path::PathBuf;

use anyhow::{anyhow, Result};
use crossterm::{
    event::{self, Event, KeyCode, KeyEvent, KeyEventKind, KeyModifiers},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use polars::prelude::{Column, DataFrame};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Constraint, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Text},
    widgets::{Block, Borders, Clear, List, ListItem, ListState, Paragraph, Wrap},
    Terminal,
};
use serde_json::{json, Value};
use tabiew::tui::{
    component::{Component, FocusState},
    table::Table as TabiewTable,
};

use crate::model::{Model, Plan};

const PAGES: [&str; 6] = [
    "Inputs",
    "Metadata",
    "Alignment",
    "Groups",
    "Settings",
    "Review",
];

struct Screen {
    model: Model,
    dry_run: bool,
    page: usize,
    focus: usize, // 0: pages, 1: actions, 2: tabiew table
    selected: usize,
    preview_quant: bool,
    table: TabiewTable,
    dialog: Option<Dialog>,
    status: String,
    dirty: bool,
    review: Option<Result<Plan, String>>,
    review_scroll: u16,
}

struct Entry {
    label: String,
    action: Action,
}

enum Action {
    MetadataFile,
    QuantFile,
    Alias(usize),
    Role(usize),
    Measure,
    QuantPreview,
    Map(String),
    Unmap(String),
    Primary,
    Secondary,
    DisplayGroup,
    Unit,
    GroupToggle(usize),
    GroupTarget(usize),
    Descriptive,
    Independent,
    Slug,
    Description,
    Seed,
    Toggle(&'static str, &'static str),
    Number(&'static str, &'static str, u32, Option<f64>),
    Method(&'static str, &'static str, Vec<&'static str>),
    Annotation(String),
    Info(&'static str),
    Refresh,
    Save,
}

#[derive(Clone)]
enum Dialog {
    Choose {
        title: String,
        choices: Vec<String>,
        selected: usize,
        purpose: ChoosePurpose,
    },
    Input {
        title: String,
        value: String,
        purpose: InputPurpose,
    },
    Quit,
}

#[derive(Clone)]
enum ChoosePurpose {
    MetadataFile,
    QuantFile,
    Alias(usize),
    Role(usize),
    Map(String),
    DisplayGroup,
    Primary,
    Secondary,
    Unit,
    Method(&'static str, &'static str),
    Annotation(String),
}

#[derive(Clone)]
enum InputPurpose {
    Verify { old: String, new: String },
    Independent,
    GroupTarget(usize),
    Slug,
    Description,
    Seed,
    Number(&'static str, &'static str, u32, Option<f64>),
}

fn frame_from_columns(names: Vec<String>, rows: Vec<Vec<String>>) -> Result<DataFrame> {
    let columns: Vec<Column> = names
        .into_iter()
        .enumerate()
        .map(|(index, name)| {
            let values: Vec<String> = rows
                .iter()
                .map(|row| row.get(index).cloned().unwrap_or_default())
                .collect();
            Column::new(name.into(), values)
        })
        .collect();
    Ok(DataFrame::new(rows.len(), columns)?)
}

fn tabiew_frame(model: &Model, page: usize, quant: bool) -> Result<DataFrame> {
    match page {
        0 if quant => {
            let suffix = format!(" Peak {}", model.measure);
            let columns = crate::data::quant_header(&model.quant_path)?;
            let rows = columns
                .into_iter()
                .map(|column| {
                    let role = if column.ends_with(&suffix) {
                        "intensity"
                    } else {
                        "feature/other"
                    };
                    vec![column, role.to_owned()]
                })
                .collect();
            frame_from_columns(vec!["Quantification header".into(), "Type".into()], rows)
        }
        2 => {
            let (metadata_only, quant_only) = model.alignment();
            let mut rows = Vec::new();
            for row in &model.table.rows {
                let old = &row[model.aliases[0]];
                let mapped = model.replacements.get(old).unwrap_or(old);
                let state = if model.replacements.contains_key(old) {
                    "VERIFIED correction"
                } else if metadata_only.contains(old) {
                    "metadata only"
                } else {
                    "matched"
                };
                rows.push(vec![old.clone(), mapped.clone(), state.to_owned()]);
            }
            for name in quant_only {
                if !model.replacements.values().any(|mapped| mapped == &name) {
                    rows.push(vec![String::new(), name, "quantification only".into()]);
                }
            }
            frame_from_columns(
                vec![
                    "Metadata filename".into(),
                    "Quant filename".into(),
                    "Alignment".into(),
                ],
                rows,
            )
        }
        3 => frame_from_columns(
            vec![
                "Source group".into(),
                "N".into(),
                "Included".into(),
                "Target".into(),
            ],
            model
                .groups
                .iter()
                .map(|group| {
                    vec![
                        group.source.clone(),
                        group.count.to_string(),
                        if group.included { "yes" } else { "no" }.into(),
                        group.target.clone(),
                    ]
                })
                .collect(),
        ),
        _ => {
            let mut names = model.table.columns.clone();
            let mut source_header = "Source group".to_owned();
            while names.contains(&source_header) {
                source_header.push('_');
            }
            names.push(source_header);
            let mut target_header = "Included / target".to_owned();
            while names.contains(&target_header) {
                target_header.push('_');
            }
            names.push(target_header);
            let rows = model
                .table
                .rows
                .iter()
                .map(|row| {
                    let source = model.source_for_row(row);
                    let group = model.groups.iter().find(|group| group.source == source);
                    let inclusion = match group {
                        Some(group) if group.included => group.target.clone(),
                        Some(_) => "excluded".into(),
                        None => "QC / blank".into(),
                    };
                    let mut out = row.clone();
                    out.push(source);
                    out.push(inclusion);
                    out
                })
                .collect();
            frame_from_columns(names, rows)
        }
    }
}

fn item(entries: &mut Vec<Entry>, label: impl Into<String>, action: Action) {
    entries.push(Entry {
        label: label.into(),
        action,
    });
}

fn value(model: &Model, section: &str, parent: &str, key: &str) -> String {
    let root = if section == "prep" {
        &model.prep
    } else {
        &model.analyses
    };
    match &root[parent][key] {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

impl Screen {
    fn new(model: Model, dry_run: bool) -> Result<Self> {
        let table = TabiewTable::new(tabiew_frame(&model, 0, false)?)
            .with_show_header(true)
            .with_extended_column()
            .striped();
        Ok(Self {
            model, dry_run, page: 0, focus: 1, selected: 0, preview_quant: false,
            table, dialog: None, status: "Tab cycles focus · Enter/Space selects · q quits · arrows navigate · table: PgUp/PgDn, ←/→".into(),
            dirty: false, review: None, review_scroll: 0,
        })
    }

    fn rebuild_table(&mut self) -> Result<()> {
        self.table = self.table.clone_with_data_frame(tabiew_frame(
            &self.model,
            self.page,
            self.preview_quant,
        )?);
        Ok(())
    }

    fn set_page(&mut self, page: usize) {
        self.page = page;
        self.selected = 0;
        self.review = None;
        if let Err(error) = self.rebuild_table() {
            self.status = error.to_string();
        }
        if page == 5 {
            self.refresh_review();
        }
    }

    fn changed(&mut self) {
        self.dirty = true;
        self.review = None;
        if let Err(error) = self.rebuild_table() {
            self.status = error.to_string();
        }
        if self.page == 5 {
            self.refresh_review();
        }
    }

    fn refresh_review(&mut self) {
        self.review = Some(self.model.plan().map_err(|error| error.to_string()));
    }

    fn entries(&self) -> Vec<Entry> {
        let mut entries = Vec::new();
        let m = &self.model;
        match self.page {
            0 => {
                if !m.is_comparison() {
                item(
                    &mut entries,
                    format!("Metadata TSV: {}", m.metadata_path.display()),
                    Action::MetadataFile,
                );
                item(
                    &mut entries,
                    format!("Quantification CSV: {}", m.quant_path.display()),
                    Action::QuantFile,
                );
                }
                if m.is_comparison() {
                    item(&mut entries, format!("Inherited metadata: {}", m.metadata_path.display()),
                        Action::Info("Metadata is fixed by the sealed study"));
                    item(&mut entries, format!("Inherited quantification: {}", m.quant_path.display()),
                        Action::Info("Quantification is fixed by the sealed study"));
                }
                item(
                    &mut entries,
                    if self.preview_quant {
                        "Show metadata table"
                    } else {
                        "Show every quantification header"
                    },
                    Action::QuantPreview,
                );
            }
            1 => {
                if !m.is_comparison() {
                for (i, role) in ["filename", "sample_id", "sample_type"].iter().enumerate() {
                    item(
                        &mut entries,
                        format!(
                            "{role} column: {}",
                            m.table
                                .columns
                                .get(m.aliases[i])
                                .map(String::as_str)
                                .unwrap_or("?")
                        ),
                        Action::Alias(i),
                    );
                }
                for (i, role) in ["sample", "QC", "blank"].iter().enumerate() {
                    item(
                        &mut entries,
                        format!("{role} role: {}", m.roles[i]),
                        Action::Role(i),
                    );
                }
                }
                if m.is_comparison() {
                    item(&mut entries, "Metadata columns and roles are fixed by the sealed study",
                        Action::Info("Inspect shared metadata in the table; corrections belong to the study"));
                }
            }
            2 => {
                if !m.is_comparison() {
                item(
                    &mut entries,
                    format!("Intensity measure: Peak {}", m.measure),
                    Action::Measure,
                );
                }
                let (missing_metadata, _) = m.alignment();
                if !m.is_comparison() {
                for name in missing_metadata
                    .iter()
                    .chain(m.replacements.keys())
                    .collect::<std::collections::BTreeSet<_>>()
                {
                    if let Some(paired) = m.replacements.get(name) {
                        item(
                            &mut entries,
                            format!("{name}  →  {paired} (remove correction)"),
                            Action::Unmap((*name).clone()),
                        );
                    } else {
                        item(
                            &mut entries,
                            format!("{name}  →  choose quant filename + VERIFY"),
                            Action::Map((*name).clone()),
                        );
                    }
                }
                }
                if m.is_comparison() {
                    item(&mut entries, format!("Peak {} and filenames inherited from study", m.measure),
                        Action::Info("Filename corrections and intensity measure belong to the study"));
                }
            }
            3 => {
                if m.is_comparison() {
                    item(&mut entries,
                        format!("Comparison primary group: {}", m.table.columns.get(m.primary).map(String::as_str).unwrap_or("?")),
                        Action::Primary);
                    item(&mut entries,
                        format!("Comparison secondary group: {}", m.secondary.and_then(|i| m.table.columns.get(i)).map(String::as_str).unwrap_or("none")),
                        Action::Secondary);
                    item(&mut entries,
                        format!("Experimental unit: {}", m.unit.and_then(|i| m.table.columns.get(i)).map(String::as_str).unwrap_or("none")),
                        Action::Unit);
                    for (i, group) in m.groups.iter().enumerate() {
                        item(&mut entries,
                            format!("[{}] {} ({} rows)", if group.included { "x" } else { " " }, group.source, group.count),
                            Action::GroupToggle(i));
                        item(&mut entries, format!("    ↳ target: {}", group.target), Action::GroupTarget(i));
                    }
                    item(&mut entries, format!("Descriptive only: {}", if m.descriptive_only { "yes" } else { "no" }), Action::Descriptive);
                    item(&mut entries, format!("Independent units: {} (type INDEPENDENT)", if m.independent { "confirmed" } else { "not confirmed" }), Action::Independent);
                } else {
                    item(&mut entries,
                        format!("Ordination display labels: {}", m.display_group.and_then(|i| m.table.columns.get(i)).map(String::as_str).unwrap_or("sample_type (default)")),
                        Action::DisplayGroup);
                }
            }
            4 => {
                item(
                    &mut entries,
                    format!("Recipe slug: {}", m.slug),
                    Action::Slug,
                );
                item(
                    &mut entries,
                    format!("Description: {}", m.description),
                    Action::Description,
                );
                item(
                    &mut entries,
                    format!("Random seed: {}", m.seed),
                    Action::Seed,
                );
                if !m.is_comparison() {
                for (parent, key, title) in [
                    ("blank_filter", "enabled", "Blank filter"),
                    ("qc_rsd_filter", "enabled", "QC RSD filter"),
                ] {
                    item(
                        &mut entries,
                        format!("{title}: {}", value(m, "prep", parent, key)),
                        Action::Toggle(parent, key),
                    );
                }
                for (parent, key, title, max) in [
                    (
                        "blank_filter",
                        "minimum_sample_to_blank_ratio",
                        "Minimum sample/blank ratio",
                        None,
                    ),
                    (
                        "blank_filter",
                        "minimum_blank_prevalence",
                        "Minimum blank prevalence (0–1)",
                        Some(1.0),
                    ),
                    (
                        "qc_rsd_filter",
                        "threshold_percent",
                        "QC RSD threshold (%)",
                        None,
                    ),
                ] {
                    item(
                        &mut entries,
                        format!("{title}: {}", value(m, "prep", parent, key)),
                        Action::Number(parent, key, 0, max),
                    );
                }
                for (parent, options) in [
                    ("missing_values", vec!["none", "half_minimum"]),
                    ("normalization", vec!["none", "total_sum", "median", "pqn"]),
                    ("transformation", vec!["none", "log2", "log10", "log1p"]),
                    ("scaling", vec!["none", "pareto", "autoscale"]),
                ] {
                    item(
                        &mut entries,
                        format!("{parent}: {}", value(m, "prep", parent, "method")),
                        Action::Method(parent, "method", options),
                    );
                }
                }
                for parent in if m.is_comparison() { &["omnibus", "differential", "plsda"][..] } else { &["pca", "pcoa"][..] } {
                    item(
                        &mut entries,
                        format!(
                            "{parent} enabled: {}",
                            value(m, "analyses", parent, "enabled")
                        ),
                        Action::Toggle(parent, "enabled"),
                    );
                }
                if !m.is_comparison() {
                for parent in ["pca", "pcoa"] {
                    for (key, role) in [("include_qc", "QC"), ("include_blanks", "blanks")] {
                        item(
                            &mut entries,
                            format!(
                                "{parent} include {role}: {}",
                                value(m, "analyses", parent, key)
                            ),
                            Action::Toggle(parent, key),
                        );
                    }
                }
                }
                for (parent, key, min) in if m.is_comparison() {
                    &[("plsda", "components", 1), ("plsda", "folds", 2), ("plsda", "repeats", 2), ("plsda", "permutations", 19)][..]
                } else {
                    &[("pca", "components", 1), ("pcoa", "components", 1)][..]
                } {
                    item(
                        &mut entries,
                        format!("{parent} {key}: {}", value(m, "analyses", parent, key)),
                        Action::Number(parent, key, *min, None),
                    );
                }
                if !m.is_comparison() {
                for (key, choices) in [
                    (
                        "stage",
                        vec!["imputed", "normalized", "transformed", "scaled"],
                    ),
                    (
                        "distance",
                        vec![
                            "bray",
                            "euclidean",
                            "maximum",
                            "manhattan",
                            "canberra",
                            "binary",
                            "minkowski",
                        ],
                    ),
                ] {
                    item(
                        &mut entries,
                        format!("PCoA {key}: {}", value(m, "analyses", "pcoa", key)),
                        Action::Method("pcoa", key, choices),
                    );
                }
                for key in m.annotation_candidates.keys() {
                    let selection = m
                        .annotations
                        .get(key)
                        .and_then(Option::as_ref)
                        .map(|p| p.display().to_string())
                        .unwrap_or_else(|| "none".into());
                    item(
                        &mut entries,
                        format!("Annotation {key}: {selection}"),
                        Action::Annotation(key.clone()),
                    );
                }
                }
            }
            _ => {
                item(
                    &mut entries,
                    "Revalidate recipe and preview exact commands",
                    Action::Refresh,
                );
                item(
                    &mut entries,
                    if self.dry_run {
                        "Finish dry-run preview (no files written)"
                    } else if self
                        .review
                        .as_ref()
                        .is_some_and(|review| review.as_ref().is_ok_and(|plan| plan.existing))
                    {
                        "Use matching hashed recipe (no files written)"
                    } else {
                        "Save new hashed recipe (no statistics run)"
                    },
                    Action::Save,
                );
            }
        }
        entries
    }

    fn chooser(
        &mut self,
        title: impl Into<String>,
        choices: Vec<String>,
        selected: usize,
        purpose: ChoosePurpose,
    ) {
        if choices.is_empty() {
            self.status = "No eligible choices available".into();
            return;
        }
        self.dialog = Some(Dialog::Choose {
            title: title.into(),
            selected: selected.min(choices.len() - 1),
            choices,
            purpose,
        });
    }

    fn input(&mut self, title: impl Into<String>, value: String, purpose: InputPurpose) {
        self.dialog = Some(Dialog::Input {
            title: title.into(),
            value,
            purpose,
        });
    }

    fn activate(&mut self, action: Action) -> Result<Option<Plan>> {
        match action {
            Action::MetadataFile => {
                let options = self.model.metadata_paths.clone();
                let current = options
                    .iter()
                    .position(|p| p == &self.model.metadata_path)
                    .unwrap_or(0);
                self.chooser(
                    "Select treated metadata TSV",
                    options.iter().map(|p| p.display().to_string()).collect(),
                    current,
                    ChoosePurpose::MetadataFile,
                );
            }
            Action::QuantFile => {
                let options = self.model.quant_paths.clone();
                let current = options
                    .iter()
                    .position(|p| p == &self.model.quant_path)
                    .unwrap_or(0);
                self.chooser(
                    "Select MZmine quantification CSV",
                    options.iter().map(|p| p.display().to_string()).collect(),
                    current,
                    ChoosePurpose::QuantFile,
                );
            }
            Action::QuantPreview => {
                self.preview_quant = !self.preview_quant;
                self.rebuild_table()?;
            }
            Action::Alias(i) => {
                let names = self.model.table.columns.clone();
                self.chooser(
                    format!(
                        "Choose {} column",
                        ["filename", "sample_id", "sample_type"][i]
                    ),
                    names,
                    self.model.aliases[i],
                    ChoosePurpose::Alias(i),
                );
            }
            Action::Role(i) => {
                let levels = self.model.role_levels();
                let selected = levels
                    .iter()
                    .position(|v| v == &self.model.roles[i])
                    .unwrap_or(0);
                self.chooser(
                    format!("Choose {} role value", ["sample", "QC", "blank"][i]),
                    levels,
                    selected,
                    ChoosePurpose::Role(i),
                );
            }
            Action::Measure => {
                let next = if self.model.measure == "height" {
                    "area"
                } else {
                    "height"
                };
                let suffix = format!(" Peak {next}");
                if !crate::data::quant_header(&self.model.quant_path)?
                    .iter()
                    .any(|column| column.ends_with(&suffix))
                {
                    return Err(anyhow!("Quantification has no Peak {next} columns"));
                }
                self.model.measure = next.into();
                self.model.replacements.clear();
                self.changed();
                self.status = format!(
                    "Peak {} selected; inspect filename alignment",
                    self.model.measure
                );
            }
            Action::Map(old) => {
                let (_, candidates) = self.model.alignment();
                self.chooser(
                    format!("Exact quantification counterpart for {old}"),
                    candidates,
                    0,
                    ChoosePurpose::Map(old),
                );
            }
            Action::Unmap(old) => {
                self.model.replacements.remove(&old);
                self.changed();
                self.status = format!(
                    "Removed correction for {old}; choose a counterpart again if necessary"
                );
            }
            Action::DisplayGroup => {
                let fields = self.model.group_fields();
                let choices = std::iter::once("sample_type (default)".to_owned())
                    .chain(fields.iter().map(|i| self.model.table.columns[*i].clone())).collect();
                let selected = self.model.display_group
                    .and_then(|i| fields.iter().position(|field| *field == i))
                    .map_or(0, |i| i + 1);
                self.chooser("Study ordination labels only; all biological rows remain", choices, selected, ChoosePurpose::DisplayGroup);
            }
            Action::Primary => {
                let fields = self.model.group_fields();
                let selected = fields
                    .iter()
                    .position(|i| *i == self.model.primary)
                    .unwrap_or(0);
                self.chooser(
                    "Primary biological grouping column",
                    fields
                        .iter()
                        .map(|i| self.model.table.columns[*i].clone())
                        .collect(),
                    selected,
                    ChoosePurpose::Primary,
                );
            }
            Action::Secondary => {
                let fields = self.model.second_fields();
                let names = std::iter::once("None".into())
                    .chain(fields.iter().map(|i| self.model.table.columns[*i].clone()))
                    .collect();
                let selected = self.model.secondary
                    .and_then(|i| fields.iter().position(|field| *field == i))
                    .map_or(0, |i| i + 1);
                self.chooser("Optional second biological grouping column", names, selected, ChoosePurpose::Secondary);
            }
            Action::Unit => {
                let fields = self.model.unit_fields();
                let names = std::iter::once("None".into())
                    .chain(fields.iter().map(|i| self.model.table.columns[*i].clone()))
                    .collect();
                let selected = self
                    .model
                    .unit
                    .and_then(|i| fields.iter().position(|field| *field == i))
                    .map_or(0, |i| i + 1);
                self.chooser(
                    "Optional experimental-unit identifier",
                    names,
                    selected,
                    ChoosePurpose::Unit,
                );
            }
            Action::GroupToggle(i) => {
                self.model.groups[i].included = !self.model.groups[i].included;
                self.model.independent = false;
                self.changed();
            }
            Action::GroupTarget(i) => self.input(
                format!("Target label for {}", self.model.groups[i].source),
                self.model.groups[i].target.clone(),
                InputPurpose::GroupTarget(i),
            ),
            Action::Descriptive => {
                self.model.descriptive_only = !self.model.descriptive_only;
                if self.model.descriptive_only {
                    for key in ["omnibus", "differential", "plsda"] {
                        self.model.analyses[key]["enabled"] = json!(false);
                    }
                }
                self.changed();
            }
            Action::Independent => self.input(
                "Type INDEPENDENT to certify independent biological units",
                String::new(),
                InputPurpose::Independent,
            ),
            Action::Slug => self.input("Recipe slug", self.model.slug.clone(), InputPurpose::Slug),
            Action::Description => self.input(
                "Recipe description",
                self.model.description.clone(),
                InputPurpose::Description,
            ),
            Action::Seed => self.input(
                "Random seed (0–2147483647)",
                self.model.seed.to_string(),
                InputPurpose::Seed,
            ),
            Action::Toggle(parent, key) => {
                let root = if ["blank_filter", "qc_rsd_filter"].contains(&parent) {
                    &mut self.model.prep
                } else {
                    &mut self.model.analyses
                };
                if self.model.descriptive_only
                    && ["omnibus", "differential", "plsda"].contains(&parent)
                {
                    return Err(anyhow!(
                        "Switch off descriptive-only mode before enabling inferential analyses"
                    ));
                }
                root[parent][key] = json!(root[parent][key] != true);
                self.changed();
            }
            Action::Number(parent, key, min, max) => {
                let section = if ["blank_filter", "qc_rsd_filter"].contains(&parent) {
                    "prep"
                } else {
                    "analyses"
                };
                self.input(
                    format!(
                        "{parent} {key} (minimum {min}{})",
                        max.map(|n| format!(", maximum {n}")).unwrap_or_default()
                    ),
                    value(&self.model, section, parent, key),
                    InputPurpose::Number(parent, key, min, max),
                );
            }
            Action::Method(parent, key, options) => {
                let section = if [
                    "missing_values",
                    "normalization",
                    "transformation",
                    "scaling",
                ]
                .contains(&parent)
                {
                    "prep"
                } else {
                    "analyses"
                };
                let current = value(&self.model, section, parent, key);
                let selected = options
                    .iter()
                    .position(|option| *option == current)
                    .unwrap_or(0);
                self.chooser(
                    format!("{parent} {key}"),
                    options.iter().map(|s| (*s).into()).collect(),
                    selected,
                    ChoosePurpose::Method(parent, key),
                );
            }
            Action::Annotation(key) => {
                let paths = self
                    .model
                    .annotation_candidates
                    .get(&key)
                    .cloned()
                    .unwrap_or_default();
                let current = self.model.annotations.get(&key).and_then(Option::as_ref);
                let selected = current
                    .and_then(|path| paths.iter().position(|p| p == path))
                    .map_or(0, |i| i + 1);
                let names = std::iter::once("None".into())
                    .chain(paths.iter().map(|p| p.display().to_string()))
                    .collect();
                self.chooser(
                    format!("Annotation {key}"),
                    names,
                    selected,
                    ChoosePurpose::Annotation(key),
                );
            }
            Action::Info(message) => self.status = message.into(),
            Action::Refresh => self.refresh_review(),
            Action::Save => {
                let plan = self.model.plan()?;
                if !self.dry_run {
                    plan.save()?;
                }
                return Ok(Some(plan));
            }
        }
        Ok(None)
    }

    fn select_choice(
        &mut self,
        purpose: ChoosePurpose,
        selected: usize,
        choices: &[String],
    ) -> Result<()> {
        match purpose {
            ChoosePurpose::MetadataFile => self.model.set_metadata(selected)?,
            ChoosePurpose::QuantFile => self.model.set_quant(selected)?,
            ChoosePurpose::Alias(i) => {
                self.model.aliases[i] = selected;
                if i == 0 {
                    self.model.replacements.clear();
                }
                self.model.refresh_groups();
            }
            ChoosePurpose::Role(i) => {
                self.model.roles[i] = choices[selected].clone();
                self.model.refresh_groups();
            }
            ChoosePurpose::Map(old) => {
                let new = choices[selected].clone();
                self.input(
                    format!("Check acquisition records; type VERIFY for {old} → {new}"),
                    String::new(),
                    InputPurpose::Verify { old, new },
                );
                return Ok(());
            }
            ChoosePurpose::DisplayGroup => {
                self.model.display_group = if selected == 0 {
                    None
                } else {
                    Some(self.model.group_fields()[selected - 1])
                };
            }
            ChoosePurpose::Primary => {
                self.model.primary = self.model.group_fields()[selected];
                if self.model.secondary == Some(self.model.primary) {
                    self.model.secondary = None;
                }
                self.model.refresh_groups();
            }
            ChoosePurpose::Secondary => {
                self.model.secondary = if selected == 0 {
                    None
                } else {
                    Some(self.model.second_fields()[selected - 1])
                };
                self.model.refresh_groups();
            }
            ChoosePurpose::Unit => {
                self.model.unit = if selected == 0 {
                    None
                } else {
                    Some(self.model.unit_fields()[selected - 1])
                };
                self.model.independent = false;
            }
            ChoosePurpose::Method(parent, key) => {
                let root = if [
                    "missing_values",
                    "normalization",
                    "transformation",
                    "scaling",
                ]
                .contains(&parent)
                {
                    &mut self.model.prep
                } else {
                    &mut self.model.analyses
                };
                root[parent][key] = json!(choices[selected]);
            }
            ChoosePurpose::Annotation(key) => {
                let path: Option<PathBuf> = if selected == 0 {
                    None
                } else {
                    self.model.annotation_candidates[&key]
                        .get(selected - 1)
                        .cloned()
                };
                self.model.annotations.insert(key, path);
            }
        }
        self.changed();
        Ok(())
    }

    fn submit_input(&mut self, purpose: InputPurpose, text: &str) -> Result<()> {
        match purpose {
            InputPurpose::Verify { old, new } => self.model.assign_filename(&old, &new, text)?,
            InputPurpose::Independent => self.model.confirm_independence(text)?,
            InputPurpose::GroupTarget(i) => {
                if text.is_empty()
                    || text == "NA"
                    || !text.as_bytes()[0].is_ascii_alphabetic()
                    || !text
                        .bytes()
                        .all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-')
                {
                    return Err(anyhow!("Target must start with a letter, use only letters/digits/_/-, and not be NA"));
                }
                self.model.groups[i].target = text.into();
                self.model.independent = false;
            }
            InputPurpose::Slug => self.model.slug = text.into(),
            InputPurpose::Description => self.model.description = text.into(),
            InputPurpose::Seed => {
                let n: u32 = text
                    .parse()
                    .map_err(|_| anyhow!("Seed must be a whole number"))?;
                if n > i32::MAX as u32 {
                    return Err(anyhow!("Seed cannot exceed {}", i32::MAX));
                }
                self.model.seed = n;
            }
            InputPurpose::Number(parent, key, min, max) => {
                let (mantissa, exponent) = text
                    .split_once(['e', 'E'])
                    .map_or((text, None), |(m, e)| (m, Some(e)));
                let valid_mantissa = !mantissa.is_empty()
                    && mantissa.bytes().any(|c| c.is_ascii_digit())
                    && mantissa.bytes().all(|c| c.is_ascii_digit() || c == b'.')
                    && mantissa.bytes().filter(|c| *c == b'.').count() <= 1;
                let valid_exponent = exponent.is_none_or(|e| {
                    let digits = e.strip_prefix(['+', '-']).unwrap_or(e);
                    !digits.is_empty() && digits.bytes().all(|c| c.is_ascii_digit())
                });
                if !valid_mantissa || !valid_exponent {
                    return Err(anyhow!("Enter a nonnegative number"));
                }
                let n: f64 = text
                    .parse()
                    .map_err(|_| anyhow!("Enter a nonnegative number"))?;
                if !n.is_finite() || n < min as f64 || max.is_some_and(|max| n > max) {
                    return Err(anyhow!(
                        "Value must be finite, at least {min}, and within its maximum"
                    ));
                }
                let section = if ["blank_filter", "qc_rsd_filter"].contains(&parent) {
                    &mut self.model.prep
                } else {
                    &mut self.model.analyses
                };
                if key == "components"
                    || key == "folds"
                    || key == "repeats"
                    || key == "permutations"
                {
                    if n.fract() != 0.0 || n > i32::MAX as f64 {
                        return Err(anyhow!("Enter a whole number <= {}", i32::MAX));
                    }
                    section[parent][key] = json!(n as u32);
                } else {
                    section[parent][key] = json!(n);
                }
            }
        }
        self.changed();
        Ok(())
    }

    fn handle_dialog(&mut self, key: KeyEvent) -> Result<Option<Plan>> {
        let Some(dialog) = self.dialog.take() else {
            return Ok(None);
        };
        match dialog {
            Dialog::Quit => match key.code {
                KeyCode::Char('y') | KeyCode::Char('Y') => return Ok(None),
                KeyCode::Char('n') | KeyCode::Char('N') | KeyCode::Esc => {
                    self.status = "Quit canceled".into()
                }
                _ => self.dialog = Some(Dialog::Quit),
            },
            Dialog::Choose {
                title,
                choices,
                mut selected,
                purpose,
            } => match key.code {
                KeyCode::Esc => self.status = "Selection canceled".into(),
                KeyCode::Up | KeyCode::BackTab => {
                    selected = selected.saturating_sub(1);
                    self.dialog = Some(Dialog::Choose {
                        title,
                        choices,
                        selected,
                        purpose,
                    });
                }
                KeyCode::Down | KeyCode::Tab => {
                    selected = (selected + 1).min(choices.len() - 1);
                    self.dialog = Some(Dialog::Choose {
                        title,
                        choices,
                        selected,
                        purpose,
                    });
                }
                KeyCode::Enter | KeyCode::Char(' ') => {
                    if let Err(error) = self.select_choice(purpose.clone(), selected, &choices) {
                        self.dialog = Some(Dialog::Choose {
                            title,
                            choices,
                            selected,
                            purpose,
                        });
                        return Err(error);
                    }
                }
                _ => {
                    self.dialog = Some(Dialog::Choose {
                        title,
                        choices,
                        selected,
                        purpose,
                    })
                }
            },
            Dialog::Input {
                title,
                mut value,
                purpose,
            } => match key.code {
                KeyCode::Esc => self.status = "Entry canceled".into(),
                KeyCode::Enter => {
                    if let Err(error) = self.submit_input(purpose.clone(), value.trim()) {
                        self.dialog = Some(Dialog::Input {
                            title,
                            value,
                            purpose,
                        });
                        return Err(error);
                    }
                }
                KeyCode::Backspace => {
                    value.pop();
                    self.dialog = Some(Dialog::Input {
                        title,
                        value,
                        purpose,
                    });
                }
                KeyCode::Char(c) if !key.modifiers.contains(KeyModifiers::CONTROL) => {
                    value.push(c);
                    self.dialog = Some(Dialog::Input {
                        title,
                        value,
                        purpose,
                    });
                }
                _ => {
                    self.dialog = Some(Dialog::Input {
                        title,
                        value,
                        purpose,
                    })
                }
            },
        }
        Ok(None)
    }

    fn handle(&mut self, key: KeyEvent) -> Result<Transition> {
        if self.dialog.is_some() {
            let quit_yes = matches!(self.dialog.as_ref(), Some(Dialog::Quit))
                && matches!(key.code, KeyCode::Char('y') | KeyCode::Char('Y'));
            self.handle_dialog(key)?;
            return Ok(if quit_yes {
                Transition::Exit(None)
            } else {
                Transition::Continue
            });
        }
        match key.code {
            KeyCode::Char('q') if key.modifiers.is_empty() => {
                if self.dirty {
                    self.dialog = Some(Dialog::Quit);
                } else {
                    return Ok(Transition::Exit(None));
                }
            }
            KeyCode::Char('[') if self.page == 5 => {
                self.review_scroll = self.review_scroll.saturating_sub(1)
            }
            KeyCode::Char(']') if self.page == 5 => {
                self.review_scroll = self.review_scroll.saturating_add(1)
            }
            KeyCode::Tab => self.focus = (self.focus + 1) % 3,
            KeyCode::BackTab => self.focus = (self.focus + 2) % 3,
            KeyCode::Up if self.focus == 0 => self.set_page(self.page.saturating_sub(1)),
            KeyCode::Down if self.focus == 0 => self.set_page((self.page + 1).min(5)),
            KeyCode::Up if self.focus == 1 => self.selected = self.selected.saturating_sub(1),
            KeyCode::Down if self.focus == 1 => {
                self.selected = (self.selected + 1).min(self.entries().len().saturating_sub(1))
            }
            KeyCode::Enter | KeyCode::Char(' ') if self.focus == 0 => self.focus = 1,
            KeyCode::Enter | KeyCode::Char(' ') if self.focus == 1 => {
                if let Some(action) = self
                    .entries()
                    .into_iter()
                    .nth(self.selected)
                    .map(|entry| entry.action)
                {
                    if let Some(plan) = self.activate(action)? {
                        return Ok(Transition::Exit(Some(plan)));
                    }
                }
            }
            _ if self.focus == 2 => {
                self.table.handle(key);
            }
            _ => {}
        }
        Ok(Transition::Continue)
    }

    fn draw(&mut self, frame: &mut ratatui::Frame<'_>) {
        let status_height = if self.page == 3 {
            9
        } else if self.page == 5 {
            5
        } else {
            4
        };
        let [body, status] =
            Layout::vertical([Constraint::Min(5), Constraint::Length(status_height)])
                .areas(frame.area());
        let [pages, menu, preview] = Layout::horizontal([
            Constraint::Length(16),
            Constraint::Percentage(40),
            Constraint::Fill(1),
        ])
        .areas(body);
        let page_items = PAGES
            .iter()
            .enumerate()
            .map(|(i, title)| {
                let style = if i == self.page {
                    Style::default()
                        .fg(Color::Black)
                        .bg(Color::Cyan)
                        .add_modifier(Modifier::BOLD)
                } else {
                    Style::default()
                };
                ListItem::new(*title).style(style)
            })
            .collect::<Vec<_>>();
        frame.render_widget(
            List::new(page_items).block(
                Block::default()
                    .title(if self.model.is_comparison() {
                        if self.focus == 0 { "Comparison pages ●" } else { "Comparison pages" }
                    } else if self.focus == 0 {
                        "Study pages ●"
                    } else {
                        "Study pages"
                    })
                    .borders(Borders::ALL),
            ),
            pages,
        );
        let entries = self.entries();
        self.selected = self.selected.min(entries.len().saturating_sub(1));
        let labels = entries
            .into_iter()
            .map(|entry| ListItem::new(entry.label))
            .collect::<Vec<_>>();
        let menu_widget = List::new(labels)
            .block(
                Block::default()
                    .title(if self.model.is_comparison() {
                        if self.focus == 1 { "Comparison decisions ●" } else { "Comparison decisions" }
                    } else if self.focus == 1 {
                        "Study decisions ●"
                    } else {
                        "Study decisions"
                    })
                    .borders(Borders::ALL),
            )
            .highlight_style(Style::default().fg(Color::Black).bg(Color::Yellow))
            .highlight_symbol("▸ ");
        let mut state = ListState::default();
        state.select(Some(self.selected));
        frame.render_stateful_widget(menu_widget, menu, &mut state);
        let (table_area, commands_area) = if self.page == 5 {
            let [table_area, commands_area] =
                Layout::vertical([Constraint::Percentage(42), Constraint::Fill(1)]).areas(preview);
            (table_area, Some(commands_area))
        } else {
            (preview, None)
        };
        let block = Block::default()
            .title(if self.focus == 2 {
                "tabiew table ●"
            } else {
                "tabiew table (all rows/columns)"
            })
            .borders(Borders::ALL);
        let inner = block.inner(table_area);
        frame.render_widget(block, table_area);
        if inner.width > 0 && inner.height > 0 {
            self.table.render(
                inner,
                frame.buffer_mut(),
                if self.focus == 2 {
                    FocusState::Focused
                } else {
                    FocusState::NotFocused
                },
            );
        }
        if self.page == 5 && self.review.is_none() {
            self.refresh_review();
        }
        let mut messages = if self.page == 5 {
            Vec::new()
        } else {
            vec![Line::from(self.status.as_str())]
        };
        if self.page == 3 {
            if self.model.is_comparison() {
                for warning in self.model.dependency_warnings() {
                    messages.push(Line::from(warning));
                }
            } else {
                messages.push(Line::from("Shared study: every biological injection is retained; grouping controls only PCA/PCoA display labels."));
            }
        }
        if self.page == 2 {
            let (meta, quant) = self.model.alignment();
            messages.push(Line::from(format!(
                "Metadata-only: {} · quantification-only: {} · corrections: {}",
                meta.len(),
                quant.len(),
                self.model.replacements.len()
            )));
        }
        if self.page == 5 {
            if let Some(review) = self.review.as_ref() {
                match review {
                    Ok(plan) => {
                        messages.push(Line::from(format!(
                            "{} {}",
                            if self.dry_run {
                                "DRY RUN"
                            } else if plan.existing {
                                "Existing matching recipe"
                            } else {
                                "Ready to save"
                            },
                            plan.folder
                                .file_name()
                                .unwrap_or_default()
                                .to_string_lossy()
                        )));
                        messages.push(Line::from(format!("SHA-256: {}", plan.hash)));
                        if let Some(area) = commands_area {
                            let text = plan
                                .commands
                                .iter()
                                .map(|command| Line::from(command.as_str()))
                                .collect::<Vec<_>>();
                            frame.render_widget(
                                Paragraph::new(Text::from(text))
                                    .wrap(Wrap { trim: false })
                                    .scroll((self.review_scroll, 0))
                                    .block(
                                        Block::default()
                                            .title(format!(
                                                "Exact commands · [/] scroll · {}",
                                                plan.folder.display()
                                            ))
                                            .borders(Borders::ALL),
                                    ),
                                area,
                            );
                        }
                    }
                    Err(error) => {
                        messages.push(Line::from(format!("Cannot save: {error}")));
                        if let Some(area) = commands_area {
                            frame.render_widget(
                                Paragraph::new(format!("Cannot save: {error}"))
                                    .wrap(Wrap { trim: false })
                                    .block(
                                        Block::default()
                                            .title("Plan validation")
                                            .borders(Borders::ALL),
                                    ),
                                area,
                            );
                        }
                    }
                }
            }
        }
        frame.render_widget(
            Paragraph::new(Text::from(messages))
                .wrap(Wrap { trim: false })
                .block(
                    Block::default()
                        .title("Status · errors")
                        .borders(Borders::ALL),
                ),
            status,
        );
        if let Some(dialog) = &self.dialog {
            draw_dialog(frame, dialog);
        }
    }
}

enum Transition {
    Continue,
    Exit(Option<Plan>),
}

fn draw_dialog(frame: &mut ratatui::Frame<'_>, dialog: &Dialog) {
    let area = frame.area();
    let width = area.width.saturating_sub(4).min(90).max(1);
    let height = match dialog {
        Dialog::Choose { choices, .. } => (choices.len() as u16 + 4)
            .min(area.height.saturating_sub(2))
            .max(3),
        _ => 7.min(area.height.saturating_sub(2)).max(3),
    };
    let box_area = Rect::new(
        area.x + (area.width - width) / 2,
        area.y + (area.height - height) / 2,
        width,
        height,
    );
    frame.render_widget(Clear, box_area);
    match dialog {
        Dialog::Choose {
            title,
            choices,
            selected,
            ..
        } => {
            let items = choices
                .iter()
                .map(|option| ListItem::new(option.as_str()))
                .collect::<Vec<_>>();
            let widget = List::new(items)
                .block(
                    Block::default()
                        .title(format!("{title} · ↑↓ Enter · Esc cancels"))
                        .borders(Borders::ALL),
                )
                .highlight_style(Style::default().fg(Color::Black).bg(Color::Yellow))
                .highlight_symbol("▸ ");
            let mut state = ListState::default();
            state.select(Some(*selected));
            frame.render_stateful_widget(widget, box_area, &mut state);
        }
        Dialog::Input { title, value, .. } => {
            frame.render_widget(
                Paragraph::new(format!("{value}▏\n\nEnter accepts · Esc cancels"))
                    .wrap(Wrap { trim: false })
                    .block(Block::default().title(title.as_str()).borders(Borders::ALL)),
                box_area,
            );
        }
        Dialog::Quit => frame.render_widget(
            Paragraph::new("Unsaved decisions will be lost. Quit? y / n")
                .wrap(Wrap { trim: false })
                .block(Block::default().title("Confirm quit").borders(Borders::ALL)),
            box_area,
        ),
    }
}

struct TerminalSession {
    terminal: Terminal<CrosstermBackend<Stdout>>,
}
impl TerminalSession {
    fn start() -> Result<Self> {
        enable_raw_mode()?;
        if let Err(error) = execute!(io::stdout(), EnterAlternateScreen) {
            disable_raw_mode()?;
            return Err(error.into());
        }
        match Terminal::new(CrosstermBackend::new(io::stdout())) {
            Ok(terminal) => Ok(Self { terminal }),
            Err(error) => {
                let _ = execute!(io::stdout(), LeaveAlternateScreen);
                let _ = disable_raw_mode();
                Err(error.into())
            }
        }
    }
}
impl Drop for TerminalSession {
    fn drop(&mut self) {
        let _ = disable_raw_mode();
        let _ = execute!(self.terminal.backend_mut(), LeaveAlternateScreen);
        let _ = self.terminal.show_cursor();
    }
}

/// Drive the in-process tabiew table alongside every menu. Writing is only possible via Review.
pub fn run(model: Model, dry_run: bool) -> Result<Option<Plan>> {
    let mut screen = Screen::new(model, dry_run)?;
    let mut session = TerminalSession::start()?;
    loop {
        session.terminal.draw(|frame| screen.draw(frame))?;
        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }
            match screen.handle(key) {
                Ok(Transition::Continue) => {}
                Ok(Transition::Exit(plan)) => return Ok(plan),
                Err(error) => screen.status = error.to_string(),
            }
        }
    }
}
