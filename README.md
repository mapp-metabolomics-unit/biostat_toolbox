# Biostats

## Install once (Ubuntu or macOS)

Use the **v0.1.0 release tag** for a reproducible toolbox checkout. Install on the machine holding the batch; a Mac user can instead SSH to an Ubuntu batch server without installing XQuartz, Chrome or R on the Mac. The Rust TUI needs Rust 1.95+, the metadata helper needs Python 3.9+ and the pinned [`mzmlwatcher`](https://github.com/oolonek/mzmlwatcher), and analyses need R 4.6.1 plus `biostats/renv.lock`.

### Ubuntu 22.04 or later

On a shared server, an **administrator** provisions the `apt` packages and system R 4.6.1 once. A researcher without sudo must **skip** those `sudo` lines: `rig list` should already show 4.6.1, and `rig run -r 4.6.1` uses that shared R from the researcher's own checkout. Do not use `rig system user-mode` on a shared server: it migrates and can remove shared R installations.

```bash
sudo apt-get update
sudo apt-get install -y curl git python3 build-essential gfortran pkg-config libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev libharfbuzz-dev libfribidi-dev
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. "$HOME/.cargo/env"
curl -LsSf https://r-lib.github.io/rig/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
sudo "$(command -v rig)" add 4.6.1
cargo install --locked --git https://github.com/oolonek/mzmlwatcher.git --rev a65011fa9489647cbdce4d9a7e84d2dcf1f2db26 mzmlwatcher
git clone --branch v0.1.0 --depth 1 https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
cd biostat_toolbox/biostats
rig run -r 4.6.1 -f bootstrap.R
cd ..
```

### macOS (when the batch is available locally)

Install the Xcode Command Line Tools if missing (`xcode-select --install`), [Homebrew](https://brew.sh/) if missing, and then:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. "$HOME/.cargo/env"
brew install r-rig
rig add 4.6.1
cargo install --locked --git https://github.com/oolonek/mzmlwatcher.git --rev a65011fa9489647cbdce4d9a7e84d2dcf1f2db26 mzmlwatcher
git clone --branch v0.1.0 --depth 1 https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
cd biostat_toolbox/biostats
rig run -r 4.6.1 -f bootstrap.R
cd ..
```

The macOS `rig add` command above uses rig's default admin mode and may request an administrator password. If a user-mode rig installation is needed, follow the [current rig user-mode instructions](https://github.com/r-lib/rig/blob/main/website/admin-vs-user-mode.qmd); older rig versions do not provide that mode.

If you already have a checkout, use the intended release tag or commit and skip `git clone`; run `bootstrap.R` from **that checkout's** `biostats/` directory. The first restore fetches pinned packages from CRAN and Bioconductor and can take a while. No root-level R environment, `params/`, web server, or browser is required.

## Run a batch

Before launching the TUI, run the helper from the toolbox root. Install [`mzmlwatcher`](https://github.com/oolonek/mzmlwatcher) on the same machine and ensure the executable is on `PATH`. By default, the helper derives the converted mzML directory as `/media/nvme8tb/msdata/<project>/<batch>` from `BATCH`, runs `mzmlwatcher scan --recursive --checksum sha256`, saves SQLite and TSV artifacts under `$BATCH/results/mzmlwatcher/`, then enriches the existing treated metadata:

```bash
BATCH="/absolute/path/to/group-repo/docs/mapp_project_XXXXX/mapp_batch_XXXXX"
python3 scripts/enrich_metadata.py --batch "$BATCH"
```

For a different mzML location, pass `--mzml-dir /absolute/path/to/converted/mzML` (or set `MZMLWATCHER_SCAN_DIR`). An explicit directory takes precedence over the environment and default. If the mounted default is unavailable, the helper also checks `$BATCH/mzml/` and mzML files directly in `$BATCH`. No prior watcher scan or TSV export is required.

Review and commit the enriched metadata in the group repository before creating recipes. The helper preserves metadata row order and existing columns, adds `injection_order` (1-based acquisition-time order), `acquisition_date`, CV-coded `polarity`, `source_raw_sha1`, `converted_mzml_sha1`, and `converted_mzml_sha256`. It requires a unique, successful watcher record with timezone-aware acquisition time and both SHA-1 values for every exact `filename` match; it refuses ambiguous timestamp ties and leaves metadata unchanged on failure. Re-running refreshes these fields; use `--output /path/to/new.tsv` instead to inspect a new metadata file without modifying the input (the TUI will not automatically select it). The watcher artifacts include local source paths; review them before publishing the batch directory.

`injection_order` is written as the first metadata column when the helper refreshes a treated TSV. Existing sealed runs do not change; a new stats run using updated metadata adds it to the saved sample table. On the dashboard's PCA, PCoA and validated PLS-DA score plots, expand **Map metadata onto sample scores** and add an `injection_order` card with **Viridis colour scale**. Add other metadata cards to explore categorical combinations without refitting the model. Expand **Legend groups** to move groups up or down and pick their colours when using categorical colouring; short legend labels retain full metadata values on hover. QC diamonds and blank crosses appear only if those injections were included in that ordination's saved recipe; PLS-DA remains biological-only.

1. **Create the shared study.** Obtain a versioned group repository (or a complete, versioned batch archive from its owner, extracted into the same directory layout). The batch must contain treated metadata, MZmine quantification and any selected annotations under `<group-repo>/docs/mapp_project_*/mapp_batch_*`. A restricted group repository requires access to reproduce its analyses; a public toolbox alone is not enough. From the toolbox root, reuse `BATCH` set above:

   ```bash
   cargo run --locked --manifest-path tools/stats-recipe-tui/Cargo.toml -- --batch "$BATCH"
   ```

   Inspect the table and verify filename discrepancies against acquisition records. Save the study recipe: its preprocessing and PCA/PCoA use **all biological samples** in the batch, regardless of later comparison targets. Its QC/blank ordination switches affect the shared axes.

2. **Validate, preview, run the study.** From the same toolbox checkout, enter `biostats/` and run the **three commands printed by the TUI**, using its actual saved paths:

   ```bash
   cd biostats
   rig run -r 4.6.1 -f run.R -- validate --dataset /path/to/saved/dataset.yaml --recipe /path/to/saved/recipe.yaml
   rig run -r 4.6.1 -f run.R -- plan --dataset /path/to/saved/dataset.yaml --recipe /path/to/saved/recipe.yaml
   rig run -r 4.6.1 -f run.R -- run --dataset /path/to/saved/dataset.yaml --recipe /path/to/saved/recipe.yaml
   ```

3. **Add comparisons when needed.** Copy the **full** `run_hash` from the study's `plan` output. Back at the toolbox root:

   ```bash
   cd ..
   STUDY_HASH=<full-64-character-study-run_hash>
   cargo run --locked --manifest-path tools/stats-recipe-tui/Cargo.toml -- --batch "$BATCH" --study "$STUDY_HASH"
   ```

   Select the biological groups and contrasts in the TUI, save, and run its three printed `validate`, `plan`, `run` commands from `biostats/`, exactly as in step 2. Repeat step 3 for additional comparisons. Each comparison is sealed separately and reuses the unchanged study preprocessing and PCA/PCoA; it does **not** rerun them. Contrasts and PLS-DA are specific to each comparison. Changing metadata, feature intensities, preprocessing or study ordination settings requires a **new study**, not enrichment of the old seal.

`run` automatically publishes `<batch>/results/stats/README.md` and the static dashboard after each sealed study or comparison. Share the **exact toolbox commit**, the batch-repository commit (or archive identifier and checksums), saved recipes, any verified metadata correction, and both study and comparison hashes. The seals record input, code, environment and output checksums; cloning a later revision can yield a different study hash. GitHub Pages links work only if the group repository publishes its `docs/` tree; deploying it may expose result files publicly.

For configuration details and scientific interpretation, see [the Biostats guide](docs/STATS_PIPELINE_BIOSTATS.md).
