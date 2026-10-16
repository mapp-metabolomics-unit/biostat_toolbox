# Batch recipe builder (tabiew-based Rust TUI)

A menu-driven terminal recipe builder embedding [tabiew](https://github.com/shshemi/tabiew)'s navigable table alongside batch-specific menus. Run it on the machine holding the batch (including over SSH); no browser, XQuartz or HTTP server is needed. Rust 1.95+ and the pinned R environment are installed using the [toolbox quickstart](../../README.md#install-once-ubuntu-server).

From the root of **your own toolbox checkout**:

```bash
BATCH="/path/to/group-repo/docs/mapp_project_XXXXX/mapp_batch_XXXXX"
cargo run --locked --manifest-path tools/stats-recipe-tui/Cargo.toml -- --batch "$BATCH"
```
This command builds a **study recipe**. The first build compiles tabiew's data engine; subsequent launches reuse the binary. Quantification is read **only through its header** and metadata appears in the interactive table. Use Inputs, Metadata, Alignment, Groups, Settings, and Review to inspect every metadata row, verify any filename correction against acquisition records (`VERIFY`), choose an optional **ordination display label** and configure shared preprocessing and PCA/PCoA. The study always includes every biological sample; its Groups page does **not** subset them.

After running the study's printed `validate`, `plan`, `run` commands from `biostats/`, copy the **full** study `run_hash` shown by `plan`. From the toolbox root, open its comparison builder:

```bash
cargo run --locked --manifest-path tools/stats-recipe-tui/Cargo.toml -- --batch "$BATCH" --study "<full-study-run_hash>"
```

This inherits the sealed study's exact inputs and corrected metadata. On Groups, choose a biological grouping column, optional second column, included groups, target labels and contrasts; on Settings choose omnibus, differential and/or PLS-DA. Review independent biological units with the study owner; confirm `INDEPENDENT` or choose descriptive-only when appropriate. PLS-DA requires at least four samples per target and folds no greater than the smallest target group. Save the comparison and run **its** three printed commands. Repeat for additional comparisons: each has its own seal and volcano/inferential outputs but reuses the study preprocessing and ordination.

**Keys:** `Tab`/`Shift+Tab` move focus between Pages, Decisions, and tabiew table; `↑`/`↓` select a page or action and `Enter`/`Space` activate it. In the table, use tabiew's arrows/`hjkl` and page keys. `Esc` closes a chooser; `q` asks before discarding changed decisions. Review shows the configuration hash, any design errors and exact commands. The TUI never runs the statistics.

The **study** Settings page offers PCA and PCoA independent **include QC** and **include blanks** switches, initially off. Selected injections participate in that shared ordination fit and appear separately from biological samples; they never enter inferential comparisons. Altering these switches or preprocessing requires a new study and new comparison recipes referencing its new hash.

`--dry-run` disables saving while retaining the interactive review and command preview. Otherwise the builder saves `<batch>/recipes/<slug>-<first-12-of-configuration-SHA-256>/` with `dataset.yaml`, `recipe.yaml` and, where required, a corrected metadata copy. The configuration fingerprint is **not** the Biostats run hash or a checksum of raw intensities. An identical configuration reuses its existing hashed directory without writing; mismatching contents are rejected. Original metadata and existing recipes are never overwritten. The R runner validates the parent study's sealed output, input, code and environment hashes before reusing it; inspect `validate` and `plan` before `run`. Each successful `run` republishes the dashboard and README under `<batch>/results/stats/`.

Generated `.yaml` files use indented JSON syntax (valid YAML), preserving categorical values such as `yes` and `no` as strings for R. Existing older recipes are not rewritten automatically.
