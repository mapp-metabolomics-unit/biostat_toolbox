# MAPP statistics refactor — paused handoff

The user requested a pause on 2026-09-29 to resume on a more powerful HPC laptop. **Do not treat this branch as released, scientifically approved, or fully verified.** No legacy files have been deleted and the historical root `renv.lock` still depends on `MAPPstructToolbox`.

## Agreed direction

- Evaluate legacy and V2 without assuming either is authoritative; retain useful V2 modules but allow reviewed numerical changes. Historical runs remain archived; new software need not read legacy run formats.
- One R analysis CLI and a read-only Shiny explorer for approved completed runs. Analysis stays offline after installation; supported targets are macOS and Linux. No Python/Rust runtime without a demonstrated need.
- Accept the usual MAPP batch tree and optional SIRIUS/CANOPUS/met-annot-enhancer annotation tables. Tolerate extra columns and peak area/height variants; require a recorded choice when both measures occur.
- Scientific identity covers input contents/mapping, methods, code, seed, R/package environment, not paths or plot styling. Runs immutable. Descriptive dashboard filters must not re-label saved inference as newly computed statistics.
- First release profile: QC, PCA/PCoA, explicit group contrasts with optional experimentally verified block, one nested-CV/permutation-validated PLS-DA, annotation candidate browsing. Unsupported designs fail explicitly. MAPP scientific owner must approve defaults and real-batch findings before replacing legacy execution.
- Internal hosted viewer limited to approved roots. No strict backward-compatible old run reader.

## Work currently in tree (uncommitted)

`git status --short` at pause reported edits to `app/app.R`, `configs/mapp_batch_00275.dataset.yaml`, `src/v2/{analysis,config,dataset,export,preprocess,runner,utils}.R`, `tests/test_v2.R`, `v2/bootstrap.R`, `v2/renv.lock`, `v2/renv/settings.json`; new `src/v2/annotations.R` and `src/v2/supervised.R`. No other files were intentionally edited. Existing user files must be preserved.

- Importer detects Peak height/area and validates metadata, unique IDs, numeric nonnegative intensity, extra columns. Recipe preflight validates methods, design rank, planned contrasts and an explicit `design.block_verified: true` for blocks.
- Runner hashes contents/mapping plus computational code and pinned environment, validates completed output checksum before reuse, stages atomically, and wires supervised + annotation modules. `mapp_v2_schema_version` is now `1.0.0`. `src/mapp_stats.R` CLI still uses `validate`, `plan`, `run`; labels and paths still say V2.
- `supervised.R` uses CRAN `pls` to implement one-hot PLS-DA with nested repeated stratified validation, fold-local preprocessing, inner-fold component selection, and label permutations. It returns `status`, `reason`, `summary`, `fold_predictions`, `permutations`, `scores`; small groups yield `withheld`. Validated is a protocol state, not evidence of significant discrimination; the UI flags near-chance/nonsignificant results.
- `annotations.R` preserves distinct candidates and evidence across three sources, including unannotated/unmatched features. Example dataset YAML now selects height and includes candidate-rich met-annot-enhancer input.
- Exporter writes TSV/RDS/plots, annotations, PLS-DA status/tables and SHA-256 `manifest$outputs`, then runner seals `COMPLETE`. The Shiny app uses `MAPP_STATS_ROOT` to list approved `results/stats_v2/<hash>` runs, verifies required object/table checksums, browses saved analyses, metadata and candidates, and offers table downloads.
- The replacement `v2/renv.lock` now includes `pls`, `ggplot2`, `yaml`, Shiny and dependencies; no V2 fork. `v2/bootstrap.R` restores pinned versions and checks the R version against the lockfile.

## Verification actually observed

- Delegate source-level/temporary smoke reports: real batch 00275 `validate`/`plan` passed (40 total samples, 2,114 features; 427 retained in its preview); temporary relocation/changed-input/ambiguous-area-height tests; synthetic PLS-DA nested validation and withholding; annotated batch normalization of 3,504 candidate records (1,568 annotated distinct features and 546 unannotated); Shiny navigation and an out-of-root symlink rejection on a sealed temporary fixture. Treat these as agent-reported, not independent release proof.
- Main ran `rig run -r 4.6-x86_64 -f bootstrap.R` from `v2/`: restore reported the library synchronized with the lockfile and completed. `renv::status()` still reports packages recorded/installed but not detected as used because source is outside the isolated `v2/` project. Resolve this warning rather than suppress it without explanation.
- Main started `rig run -r 4.6-x86_64 -f ../tests/test_v2.R` but **cancelled it at the user's request**; no result. A main `rig run ... -f run.R validate --dataset ...` failed at `rig` argument parsing before the CLI ran: `rig` requires `--` before dash-prefixed forwarded arguments. Both jobs are stopped. No main end-to-end run, macOS UI visual verification, Linux installation, or scientific approval has been completed.

## Known integration risk — fix first on resumption

`src/v2/utils.R:write_tsv()` was changed to correctly quote fields, but `app/app.R:saved_table()` currently calls `read.delim(..., quote = "")`. The app will misread exported quoted headers/values. Change the reader to handle `write.table`'s quoted TSV (standard `quote = "\""`), then rerun a real completed-run UI smoke. This was noticed immediately before the stop request and deliberately left unedited.

Other checks: inspect renamed/removed-file implications of output SHA sealing, whether metadata-only samples should be excluded rather than hard-failed for representative batches, PCA/PCoA one-dimensional static exports, and whether the installed `pls` version is available on Linux. PLS-DA performance and scientific defaults require owner review. The root legacy runner, old consumers (`plot_selected_boxplots.R`, `analyze_spectral_modules.R`, `reprocess_stats.R`), README and installation instructions have **not** been cut over; root fork dependency remains for historical use.

## Concrete resumption sequence

1. Fix the TSV reader mismatch and review the affected new source files. From repo root: `git status --short`; inspect `src/v2` and `app/app.R` before edits. Preserve uncommitted changes.
2. Run `cd v2 && rig run -r 4.6-x86_64 -f ../tests/test_v2.R` (macOS version alias observed here); repair failures. On another machine, use its actual installed R 4.6.1 executable or alias. The test asserts area/height ambiguity, portable hashes, core contrasts and immutable-run corruption detection.
3. Run actual batch 00275 through CLI `validate`, `plan`, then `run` using the corrected `rig` forwarding syntax (likely `rig run -r 4.6-x86_64 -f run.R validate -- --dataset ../configs/mapp_batch_00275.dataset.yaml --recipe ../configs/mapp_batch_00275.recipe.yaml`; verify syntax with `rig run --help`) or invoke R 4.6.1 `Rscript run.R validate --dataset ...` directly from `v2/`. Do not claim this passed yet. Evaluate PLS-DA separately with reasonable compute/permutation settings, not on the current weak laptop.
4. Launch Shiny with `MAPP_STATS_ROOT=<approved docs root>` under the V2 R library and navigate the real run in a browser: catalogue, PCA/PCoA, contrasts, PLS-DA validation/withheld state, candidate evidence, metadata filters and verified downloads. Verify privacy/symlink behavior.
5. Exercise at least one other representative real batch, including role/annotation differences, and synthetic dual-measure tables. Resolve installation lock/status and verify Linux on the HPC system.
6. Only after smoke evidence: document supported CLI/app workflow and scientific limitations; inventory real consumers; remove or archive superseded code and fork from the supported install path. Do **not** retire legacy execution or claim scientific equivalence without MAPP scientific-owner signoff.

This file exists only because the user explicitly requested a durable handoff. No cleanup/deletion is authorized merely by these notes.
