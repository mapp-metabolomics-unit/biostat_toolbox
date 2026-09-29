# MAPP statistical pipeline V2

## Status

V2 is a side-by-side replacement under active development. The legacy entry point, `src/biostat_toolbox.r`, and its root `renv.lock` remain authoritative for reproducing historical runs. V2 never imports `MAPPstructToolbox`.

The first vertical slice implements validated MZmine import, blank filtering, optional QC-RSD filtering, imputation, normalization, transformation, scaling, PCA, PCoA, omnibus linear-model tests, planned contrasts, static exports, content-addressed runs, and a generic read-only Shiny explorer.

## Commands

From the repository root using the existing environment:

```bash
Rscript src/mapp_stats.R validate --dataset configs/mapp_batch_00275.dataset.yaml --recipe configs/mapp_batch_00275.recipe.yaml
Rscript src/mapp_stats.R plan --dataset configs/mapp_batch_00275.dataset.yaml --recipe configs/mapp_batch_00275.recipe.yaml
Rscript src/mapp_stats.R run --dataset configs/mapp_batch_00275.dataset.yaml --recipe configs/mapp_batch_00275.recipe.yaml
```

`validate` checks matrix/metadata alignment, roles, group replication, and method prerequisites. `plan` prints hashes, the output path, and the number of features each preprocessing filter would remove without changing the batch. Review this preview before `run`, especially when QC-RSD or blank filtering removes a large fraction of features. `run` writes an immutable completed run under `results/stats_v2/<run_hash>`.

To open a completed run in the generic explorer, use the isolated V2 environment:

```bash
cd v2
MAPP_STATS_RUN=/absolute/path/to/results/stats_v2/<run_hash> Rscript app.R
```

## Run identity

V2 records three SHA-256 identities:

- `dataset_hash`: dataset ID plus checksums of metadata, quantification, and annotation inputs.
- `preprocess_hash`: dataset hash, roles, preprocessing recipe, V2 source checksums, R/package versions, root lockfile checksum, and Git state.
- `run_hash`: preprocessing hash, design, analyses, seed, and environment identity.

Output paths and presentation settings are excluded from scientific hashes. Vector order is preserved because it may encode contrast or model order. Runs are assembled in a staging directory, published atomically, protected by a per-hash lock, and marked with `COMPLETE` only after all exports succeed.

## Processing order and matrix branches

1. Import raw peak heights while retaining samples, QCs, and blanks.
2. Apply blank filtering before removing any sample role.
3. Restrict the inferential branch to biological samples.
4. Impute nonpositive/missing values.
5. Normalize between samples.
6. Optionally evaluate normalized QC RSD and remove unstable features.
7. Transform.
8. Scale.

Analyses deliberately use different branches:

- PCA consumes `scaled`.
- Bray-Curtis PCoA consumes nonnegative `normalized` values.
- Linear models and volcano plots consume `transformed`, not Pareto-scaled, values.
- Feature exploration can display any saved stage.

Every volcano plot is rendered from `tables/differential.tsv`; the interactive and static coordinates therefore cannot diverge.

## Statistical interpretation

The omnibus model tests whether the group factor explains variation in each feature. Planned contrasts report an effect, its standard error, test statistic, raw p-value, and Benjamini-Hochberg q-value. With a `log2` transformation, the effect is a log2 fold change.

A block can be declared as `design.block`. It must only be used when its levels represent genuine matching, pairing, plate, batch, or another experimental blocking structure. V2 does not infer this silently.

PLS-DA is intentionally not part of the first slice. It will be added through official `structToolbox` only with repeated stratified cross-validation and permutation testing; an unvalidated scores plot will not be presented as evidence.

## Batch 00275 decisions

The quantitative matrix contains all 30 samples, seven QCs, and three blanks. Because no blank-subtraction module is recorded in its MZmine batch, V2 blank filtering is enabled. The initial thresholds are visible recipe choices, not hidden defaults.

`ATTRIBUTE_mutant` has six balanced levels with five samples each. The omnibus group test is supported. Only `WT_O_vs_WT_N` is initially declared. Mutant-versus-control contrasts should be added after confirming which wild type is the appropriate reference.

The `ATTRIBUTE_replicate` values A-E occur across every group. The recipe does not block on replicate until the experimental owner confirms that these are matched blocks.

## Dependency policy

The root environment retains the pinned `MAPPstructToolbox` fork solely for legacy runs. V2 uses a separate environment and the official Bioconductor `structToolbox`, plus focused maintained packages when their implementation is clearer or more robust. Namespace-qualified calls prevent class collisions.

Current official Bioconductor releases require a newer R than the legacy R 4.2.2 environment. Install R 4.6 side by side with `rig add 4.6`; then enter `v2/` and use `rig run -r 4.6 -f bootstrap.R` to create and snapshot the V2 lockfile. Continue using `rig run -r 4.6` for V2 commands so the default legacy R is unchanged. Do not update the root lockfile as part of that operation.

## Next implementation stages

1. Add QC drift correction after `injection_order` validation.
2. Add official `structToolbox` PLS-DA with nested/repeated validation and permutation tests.
3. Add preprocessing previews to Shiny without committing a run.
4. Add background job submission; submitted jobs must execute the same CLI recipe.
5. Add annotation/NPC joins and composition views to the generic explorer.
6. Regression-test batches 00196, 00270, and 00275 before retiring any legacy output.
