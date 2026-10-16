mod data;
mod model;
mod ui;
mod study;

use std::env;
use std::io::{self, IsTerminal};
use std::path::PathBuf;

use anyhow::{bail, Result};

fn run() -> Result<()> {
    let mut args = env::args_os().skip(1);
    let mut batch = None;
    let mut study = None;
    let mut dry_run = false;
    while let Some(arg) = args.next() {
        if arg == "--batch" {
            if batch.is_some() {
                bail!("--batch may be supplied only once");
            }
            batch = Some(PathBuf::from(
                args.next()
                    .ok_or_else(|| anyhow::anyhow!("--batch requires a path"))?,
            ));
        } else if arg == "--study" {
            if study.is_some() { bail!("--study may be supplied only once"); }
            study = Some(args.next().ok_or_else(|| anyhow::anyhow!("--study requires a full study hash"))?
                .into_string().map_err(|_| anyhow::anyhow!("--study hash must be UTF-8"))?);
        } else if arg == "--dry-run" {
            dry_run = true;
        } else if arg == "--help" || arg == "-h" {
            println!("Usage: stats-recipe-tui --batch PATH [--study FULL_STUDY_HASH] [--dry-run]\nWithout --study, create a shared batch-wide study. With --study, create an incremental comparison from a sealed study. No statistics are run.");
            return Ok(());
        } else {
            bail!(
                "Unrecognized argument {:?}; expected --batch PATH [--study FULL_STUDY_HASH] [--dry-run]",
                arg
            );
        }
    }
    let batch = batch.ok_or_else(|| anyhow::anyhow!("--batch PATH is required"))?;
    if !io::stdin().is_terminal() || !io::stdout().is_terminal() {
        bail!("An interactive terminal is required. Run this on the batch host over SSH, not through a pipe.");
    }
    let model = if let Some(hash) = study {
        model::Model::load_comparison(&batch, &hash)?
    } else {
        model::Model::load(&batch)?
    };
    if let Some(plan) = ui::run(model, dry_run)? {
        println!("\nRecipe configuration SHA-256: {}", plan.hash);
        println!("Recipe directory: {}", plan.folder.display());
        println!("\nValidate and plan first; run only after inspecting their output:");
        for command in &plan.commands {
            println!("{command}");
        }
        if dry_run {
            println!("\nPreview only: no files written and no statistics run.");
        } else if plan.existing {
            println!("\nUsing the matching existing recipe; no files written or statistics run.");
        } else {
            println!(
                "\nSaved a new hashed recipe. Original metadata unchanged. No statistics run."
            );
        }
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("Stopped: {error:#}");
        std::process::exit(1);
    }
}
