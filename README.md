# biostat_toolbox

> The historical workflow remains available through `src/biostat_toolbox.r`. A new, side-by-side pipeline is documented in [docs/STATS_PIPELINE_V2.md](docs/STATS_PIPELINE_V2.md). V2 uses immutable content-addressed runs and does not replace the legacy runner or its environment.

Script-based metabolomics analysis utilities for the MAPP biostat workflow.

## Browse completed V2 runs without Shiny

Batch 00275's previous ten V2 runs and copied derived inputs were cleared on request. New owner-produced full-batch and merged-group runs are sealed and published together in one reusable static dashboard under `results/stats_v2/dashboard/`; the catalogue currently indexes six completed runs. The viewer selects sealed runs by hash from `runs.json`; publishing another run refreshes that catalogue rather than making a separate dashboard. To refresh and preview it:

```bash
cd v2
rig run -r 4.6.1 -f ../scripts/publish_stats_dashboard.R -- --batch ../../didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275
cd ..
python3 -m http.server 8765 --bind 127.0.0.1 --directory ../didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/results/stats_v2
```

Open `http://127.0.0.1:8765/dashboard/` after a completed CLI run and publication. Use the run selector or bookmark `?run=<64-character-run-hash>`. Rerun the publishing command after each completed CLI run to refresh the index. The publisher requires the pinned V2 R packages, but **readers need only an HTTP(S) browser**. The Python server above is for local preview; the dashboard does not work from `file://` and requires HTTPS (or loopback HTTP) for browser SHA-256 verification. Static hosting has no built-in login: protect the dashboard **and every run file** behind the same access controls. The browser verifies saved TSV/PNG/PDF checksums and does not recompute statistics. Structure images use the [Natural Products API](https://api.naturalproducts.net/latest/docs#/depict/depict_2d_molecule_depict_2D_get): hovering a feature point or viewing its gallery sends saved SMILES and the reader's network address to that external service; when depiction fails, saved SMILES remains available.

**Provenance & downloads → Completed jobs** lists sealed V2 analyses with timestamps, short descriptions, a stable recipe SHA-256 and the exact input, mapping, code or environment fields changed since the preceding run for that dataset. Same recipe fingerprint does **not** imply identical source data or implementation; those differences intentionally change the run hash. Filter, sort or download the list, or select a run hash to inspect it. Publishing rebuilds `results/stats_v2/params_log.tsv` from sealed manifests; the historical `results/stats/params_log.tsv` is separate and was not cleared.

**Advanced · new analysis** is available once a canonical full-batch run is published. It starts from that run's checksum-verified complete metadata, even when the picker shows a subset. Choose a categorical biological grouping field; include, exclude or merge original groups, preview counts, configure V2 filters and analyses, give the job a short description, then download a JSON plan. The local materializer verifies the sealed source and canonical input SHA-256 and writes only dataset/recipe YAML under the **user's batch** `recipes/<slug>/`; it never copies metadata or quantification. The previous `stats_v2_inputs/` derived-input archive and prior V2 runs have been cleared as requested. The browser does not execute statistics or alter a saved run. See [V2 recipe instructions](v2/README.md#build-a-new-recipe-from-the-static-dashboard).

**Overview → Sample metadata** displays every saved sample and every field from the selected run’s checksum-verified `sample_metadata.tsv`. Filter or sort columns, choose visible columns and download the displayed TSV; these are read-only views of the saved rows.

**Overview → Feature retention & annotation** reconciles saved raw feature IDs with checksum-verified `blank_filter.tsv`, `qc_rsd.tsv` and retained features. It shows the count removed by each sequential filter, saved thresholds, and a filterable/downloadable audit for every raw feature. Click an ID to inspect its raw intensities. Annotation candidate records are not unique features and do not decide filtering.

The **PLS-DA** tab displays checksum-verified held-out accuracy, label-permutation evidence and fold predictions when a completed run includes PLS-DA. Validated runs also display descriptive full-fit sample scores, feature loadings and VIP values; these are **not** independently validated feature importance. An OPLS S-plot is not produced by this PLS-DA workflow. Batch 00275's published merged-group run reports held-out balanced accuracy 0.64 and a 99-permutation p-value of 0.01; these saved validation results do not establish chemical identities or external generalization.

The **Feature explorer** browse table shows one row per feature; source-specific candidate records remain selectable within that feature's evidence panel. **Saved InChIKey2D agreement** and **Saved agreeing structure sources** are taken directly from `annotation_horizontal.tsv` (`sources_number_IK2D` and `sources_IK2D`); the viewer does not calculate consensus from candidate rows or SMILES. GNPS, ISDB and SIRIUS can contribute structure evidence; CANOPUS reports class annotations, not structure consensus. **Candidate record source** matches features with saved SIRIUS, CANOPUS or `met_annot_enhancer` candidate rows; it does not recalculate the saved agreement. GNPS can contribute to agreement without a candidate record. Reported-taxon and text searches ignore letter case.

The **Chemical classes** tab drills through checksum-verified CANOPUS NPC pathway → superclass → class → measured feature. Choose any available categorical biological metadata field and two groups to compare; a filterable/downloadable table and paired bars show each branch’s mean *raw* intensity per sample, B/A ratio and percentage of its parent’s annotated raw intensity. Feature IDs open the inline Feature explorer. The comparison is descriptive: CANOPUS classes are predictions, unannotated features are excluded, and no new statistical test or normalization is performed.

Hover source choices and the depiction gallery treat `met_annot_enhancer` SMILES as ISDB structures, not a fourth structural source; the candidate-record provenance remains `met_annot_enhancer`. Depictions are offered only for saved GNPS, ISDB or SIRIUS structures.

Every dashboard table has case-insensitive per-column text filters, ascending/descending column sorting, and a **Columns** chooser. Numeric columns also accept `>2`, `>=2`, `<2`, `<=2`, `=2`, or an inclusive range such as `1..3` in their filter box; for example, filter **Saved InChIKey2D agreement** with `>2` to show values above two. **Download displayed TSV** exports all rows matching the current filters in the chosen order and visible columns, including rows beyond the current page; it is a descriptive view, not a replacement for checksum-verified saved TSV files. Full HTTP(S) URLs in saved evidence, including `organism_wikidata` and `identification_links`, are clickable and open in a new tab.

The **Volcano** tab shows one card per group. Drag cards into **Left · reference** and **Right · numerator**, or select a card and activate a target with the keyboard; **Swap sides** reverses the view. Only saved group comparisons can be plotted. Reversing a saved comparison negates its displayed effects while retaining saved p/q-values; downloadable original PNG/PDF plots retain their saved orientation. The x-axis is symmetric around zero so the zero fold-change line stays centered for every comparison.

Hovering a feature point on **Volcano**, **PCA or PLS-DA feature loadings**, or a **Feature explorer intensity plot** shows its saved feature ID, m/z and retention time, and requests a 2D depiction when a saved SMILES is within the service's 5000-byte limit. In the left **Settings · feature hover** menu, choose which saved metadata fields appear and which source supplies the depiction. Selecting a source never falls back to a different source for a feature lacking its SMILES. Hover depictions use the same external API as the feature gallery; hovering sends that saved SMILES and your network address to the service.

Clicking a feature point in **PCA or PLS-DA feature loadings**, **Volcano**, or the **PLS-DA VIP** plot opens the same Feature explorer structures, measured metadata, candidate records and descriptive intensity plot in a foldable side panel, leaving the source plot visible. Collapse the panel to restore plot width, or open the full Feature explorer with the selected feature retained.
Select **Detach to window** in the feature side panel to move its detail view onto another screen. The dashboard stays on the current plot; feature selections, run changes and group colours update the separate window, and choosing a feature there updates the main explorer. Use **Reattach here** in the dashboard or **Reattach to dashboard** in the separate window to restore the inline panel. Allow popups for the dashboard if the browser blocks the new window; the inline panel remains available when blocked. Both windows read the same published, checksum-verified run files from the same origin.

**PCA and PLS-DA feature loadings** default to teal for features with a saved GNPS, ISDB or SIRIUS SMILES (including enhancer-reported ISDB), and gray otherwise. In **Display options · feature colours**, both loading views and Volcano map categorical fields such as CANOPUS NPC pathway to editable colour groups. Numeric fields use a continuous Viridis colourbar instead: saved 2D-key agreement (`sources_number_IK2D`, 0–3), saved ISDB taxonomic score (`isdb_score_taxo`, 0–8), and numeric measured feature fields such as m/z and retention time. Zero is recorded, not missing; **Not recorded** and **Multiple values** remain separate, editable exception colours. **Candidate record source** describes candidate rows, not the number of agreeing structures. Volcano still uses diamonds for saved q ≤ 0.05.

In 3D, rotating the PCA sample-score plot rotates its feature-loading plot to the same camera and vice versa. Bray-Curtis PCoA has saved sample coordinates but no native PCA-style feature loadings; the distance ordination does not define feature coefficients.

Sample grouping and role controls are local to **Descriptive sample intensities**, not the sidebar. The sidebar's single **Group colours · all plots** picker sets a consistent colour for each sample group in PCA, PCoA, PLS-DA and descriptive feature intensities. In the intensity card, open **Display options** to choose the saved intensity stage, groups and sample roles, and facet by a categorical sample-metadata field. The plot stays visible when the controls are collapsed. Stage, grouping, role, facet, selected groups and menu state are saved in this browser per sealed run: they survive changing features, moving between the inline drawer and full explorer, reloading the page and returning to that run. A detached same-origin window shares these settings. An unseen run starts with its defaults; clearing site storage clears saved display preferences. Facets share the intensity scale; these choices change only the descriptive plot, not saved statistics.

A public GitHub Pages address cannot fetch results from a private GitHub repository on behalf of its visitors. Even a Pages site *built from a private source repository* can be public, exposing every published data asset. [Private GitHub Pages access control](https://docs.github.com/en/enterprise-cloud@latest/pages/getting-started-with-github-pages/changing-the-visibility-of-your-github-pages-site) requires an eligible GitHub Enterprise Cloud organization and protects the published site only when configured as private. Otherwise serve the viewer and data together behind authentication on your own host; do not put private run files in a public Pages deployment.

## View a completed V2 run

The R 4.2 installation instructions below are for **historical runs**, not the V2 viewer. For the example batch 00275, from this repository's root (with its `didier-reinhardt-group` sibling checkout and R 4.6.1 installed):

```bash
cd v2
rig run -r 4.6.1 -f bootstrap.R  # first time only
MAPP_STATS_ROOT=../../didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275 \
  rig run -r 4.6.1 -e 'shiny::runApp("../app", host="127.0.0.1", launch.browser=FALSE)'
```

Keep the last command running and open the **Listening on http://127.0.0.1:PORT** URL printed in its terminal. Shiny chooses a free port; do not assume port 43117. If running on `commons-server` over SSH, [forward that port from your laptop](docs/STATS_PIPELINE_V2.md#access-from-another-machine)—the server's `127.0.0.1` is not your laptop's.
`MAPP_STATS_ROOT` points to the **batch directory**, not a run hash; substitute your own batch directory for other inputs. The newest completed run is selected on opening. Explore PCA, PCoA, **Volcano plots** and **Feature explorer**. Click a PCA *loading* point or a volcano point to open the right-side collapsible Feature explorer drawer; feature IDs in the explorer's table and picker instead select a feature in the full-width page. The gallery depicts one representative SMILES per saved matching 2D InChIKey across GNPS, ISDB and SIRIUS; differing or unknown keys remain separate, and expandable details retain each source's original SMILES. Named stars indicate which of the three sources share that key, not confirmed chemical identity. An ISDB star fades from pale yellow to bright green as its saved taxonomic score increases on the displayed 0–8 visual scale, with score and reported taxon in its tooltip; a missing score leaves it neutral. A descriptive intensity box plot sits beside the gallery by default; all saved candidate records for that feature stay directly below, even when browse filters exclude them. Met-annot-enhancer candidates show reported taxa and saved Wikidata links. The browser sends depicted SMILES to the public Natural Products depiction API. The viewer never changes immutable statistical results. Saved PNG/PDF plots and TSV tables remain under `<batch>/results/stats_v2/<run_hash>/`; see the [V2 viewer walkthrough](docs/STATS_PIPELINE_V2.md#read-only-explorer) for filters and downloads.

The dataset and metadata sidebar is hidden in **Feature explorer** to give structures, intensity and matching candidate rows the full width; switch to another tab to adjust its saved intensity stage, grouping or sample-role controls. **Feature explorer** retains all feature and annotation filters, including source, consensus, NPC classifications and exact or partial-name **Reported taxon (ISDB reponderated only)** selection. Search the compact scrollable candidate table below the selected feature, page through every saved candidate and choose its visible columns; clicking a candidate row focuses that source record without clearing filters. Reported taxa are annotation evidence, not confirmed sample origins; taxon filters never filter samples or the selected feature's full structures and pinned records.

Measured features removed by preprocessing retain saved raw intensities: the default-on, optional intensity plot and plot-driven drawer show raw peak height/area when the selected stage lacks a feature and explain the saved filter decision. The newest run also includes input measurement metadata. Batch 00275 contains peak **height** only.

The batch-owned full recipe omits `design.contrasts`, which defaults to all 15 pairwise volcanoes for its six original groups without merging them. Provide an explicit contrast list to restrict them; a rerun with identical effective recipe, source checksums and code/environment reuses the same completed run hash.

The sections below describe the historical workflow.

## Supported platforms

`biostat_toolbox` is supported on `macOS` and `Linux`.
Windows is not supported.

## Installation

### System prerequisites

- `R 4.2.x`
- A working compiler toolchain (Xcode CLT on macOS, `build-essential` on Linux)
- `pandoc` on `PATH`

### Steps

```bash
git clone https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
cd biostat_toolbox
Rscript install.R
cp params/params_template.yaml params/params.yaml
cp params/params_user_template.yaml params/params_user.yaml
```

Edit the two param files for your dataset, then verify:

```bash
Rscript scripts/check_install.R
```

### Run

```bash
Rscript src/biostat_toolbox.r
```

You can also run the main script with explicit parameter files:

```bash
Rscript src/biostat_toolbox.r \
  --params /path/to/params.yaml \
  --params-user /path/to/params_user.yaml
```

PCA, PLSDA, and PCoA share the `ordination` plot styling and sample-label settings from `params.yaml`.

`npc_summed_intensity` can export a boxplot comparing the summed raw intensity of all retained features assigned to selected CANOPUS NPC terms. Leave all term lists empty to skip the plot.

```yaml
npc_summed_intensity:
  pathway: []
  superclass:
    - "Carotenoids (C40)"
  class:
    - "Triacylglycerols"
    - "Diacylglycerols"
  min_probability: 0
  transform: "log10"
  raw_export: TRUE
  ratios:
    enabled: TRUE
    denominator_level: "pathway"
    pseudocount: 0
```

When enabled, compact outputs are organized under `NPC_summed_intensity/filtered` and written from `DE` after filters and scaling. If `raw_export` is `TRUE`, matching outputs are written under `NPC_summed_intensity/raw` from `DE_original` intensities and features but using the filtered sample set, so blanks/QCs excluded by the params stay excluded. Combined PDF/PNG figures and complete intensity, statistics, ratio, and feature-driver TSV files are retained. Per-term dashboards and the legacy `NPC_feature_explorer.html` are no longer generated; interactive exploration is consolidated in `data_explorer.html`. The denominator pathway for ratio exports is resolved from the NP-Classifier taxonomy dictionary, not inferred from the observed CANOPUS pathway column.

## Reprocess Existing Stats Runs

Use `src/reprocess_stats.R` to rerun all `biostat_toolbox` stats outputs from archived result folders. Each source folder must contain the original `params.yaml`; if `params_user.yaml` is present, it is reused as the base configuration. Pass `--params-user` after moving a project to recursively override archived user settings such as `paths.docs` and `operating_system.pandoc`. The command-line output root always replaces `paths.output`.

```bash
Rscript src/reprocess_stats.R \
  --stats-dir /path/to/results/stats \
  --output-root /path/to/new/results/stats \
  --params-user params/params_user.yaml \
  --dry-run
```

Remove `--dry-run` to launch the reprocessing. Outputs are written under params-derived new hash folders in `--output-root`, and `reprocess_manifest.tsv` maps each original hash to the new hash. By default, reprocessing imports missing plot-default sections (`ordination`, `npc_summed_intensity`) from `params/params.yaml`; use `--defaults-yaml /path/to/params.yaml` to choose another defaults file.

To apply shared setting updates, copy and edit `params/params_override_template.yaml`, then pass it as an override YAML:

```yaml
by_target:
  ATTRIBUTE_part:
    colors:
      all:
        key:
          - Yellowpatch
          - Greenpatch
          - Uppwing
          - Lowwing
          - Body
        value:
          - "#C99700"
          - "#2E8B57"
          - "#E76F51"
          - "#4D96D7"
          - "#7B4FA3"
  ATTRIBUTE_sex:
    colors:
      all:
        key:
          - female
          - male
        value:
          - "#B55AA0"
          - "#3E78B2"
```

Color overrides are treated as a palette: only keys already present in a source `params.yaml` are updated. This keeps mixed batches safe when different archived runs compare different group subsets.

```bash
Rscript src/reprocess_stats.R \
  --stats-dir /path/to/results/stats \
  --output-root /path/to/new/results/stats \
  --override-yaml params/params_override.yaml
```

Useful options:

- `--include` / `--exclude`: comma-separated original hashes to select or skip
- `--params-user`: merge current machine paths and operating-system settings over archived user parameters
- `--overwrite`: rerun even if the predicted output already contains `session_info.txt`
- `--stop-on-error`: stop at the first failed run

## Data Explorer

Use `src/generate_data_explorer.R` for interactive feature and NPC exploration. It provides feature, molecular-component, NPC class, superclass, and pathway views; sample metadata filtering; grouping, faceting, and coloring; raw, log-transformed, and sample-relative intensities; annotation-score filters; side-by-side composition treemaps; archived-statistics and dynamic Welch-test volcano modes; and nested drill-downs to individual features. Archived Volcano mode imports completed contrasts from `reprocess_by_extraction_phase`, matches them by grouping variable, groups, and sample set, and displays their stored fold-change and p-value coordinates exactly. When no exact archived sample set is available, the panel explicitly identifies its Dynamic Welch fallback.

```bash
Rscript --vanilla src/generate_data_explorer.R \
  --params /path/to/params.yaml \
  --params-user /path/to/params_user.yaml \
  --output-dir /path/to/output
```

The explorer writes `data_explorer.html` and its lazy-loaded assets to `data_explorer_assets`. Open it through a local web server because browsers commonly block chunk loading from a direct `file://` URL. Inferential NPC statistics and pathway-normalized ratios remain available in the compact TSV exports from the main workflow.

## Selected Boxplots

Use `src/plot_selected_boxplots.R` to export individual plots for a chosen set of features from an existing `biostat_toolbox` result set.

### Basic usage

```bash
Rscript src/plot_selected_boxplots.R \
  --hash 1c63fd85afc717643ca8186ff99e6ae4 \
  --features "123,456,789"
```

This command:

- reads `params/params.yaml` and `params/params_user.yaml` by default
- uses the configured stats directory from the param files and appends the hash passed with `--hash`
- writes plot images into the resolved hash-specific result directory by default
- writes the exported intensity table to `selected_boxplot_data.csv` inside the output directory

### Common options

```bash
Rscript src/plot_selected_boxplots.R \
  --hash 1c63fd85afc717643ca8186ff99e6ae4 \
  --features "123,456,789" \
  --plot-type violin_box \
  --output-dir /tmp/selected_boxplots
```

- `--hash`: hash of the result directory to read under the configured stats directory
- `--features`: comma-separated feature IDs to plot
- `--features-file`: text file with one feature ID per line
- `--results-dir`: explicitly point to a result directory containing `DE.rds` and `foldchange_pvalues.csv`
- `--output-dir`: override the default plot output directory
- `--plot-type`: one of `box`, `violin`, or `violin_box`
- `--data-output`: override the path of the exported long-format table
- `--pvalue-column`: choose the p-value column when `foldchange_pvalues.csv` contains more than one

### Notes

- Run `Rscript src/biostat_toolbox.r` first so the expected result files exist.
- If multiple result hashes exist for the same dataset, use `--hash` to choose the one to plot from.
- The default output location is the resolved result directory, which comes from `params_user.yaml` and the selected hash.
- If you need non-default parameter files, use `--params` and `--params-user`.

## CANOPUS NPC Heatmap

Use `src/canopus_npc_heatmap.R` to build a standalone interactive feature heatmap from a SIRIUS `canopus_structure_summary.tsv`, the paired MZmine quant table, and the treated sample metadata.

The script uses the repo's existing R stack:

- `iheatmapr` for the interactive heatmap and annotation tracks
- `plotly/htmlwidgets` through `iheatmapr` for browser-based interaction and HTML export
- `readr/dplyr/tidyr/stringr` for the data joins and preprocessing

### Basic usage

```bash
Rscript src/canopus_npc_heatmap.R \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv
```

By default this command:

- infers the MZmine quant table as `results/mzmine/<batch>_quant.csv`
- infers the metadata table as `metadata/treated/<batch>_metadata.tsv`
- keeps only `sample_type == sample`
- selects the top 300 CANOPUS-annotated features by summed intensity
- applies `log10(x + 1)` followed by row z-scoring
- writes `canopus_npc_feature_heatmap.html` next to the CANOPUS file
- writes a companion TSV with the exact plotted data matrix and annotations

### Example with explicit outputs

```bash
Rscript src/canopus_npc_heatmap.R \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --output-html /tmp/canopus_npc_feature_heatmap.html \
  --top-n 500 \
  --cluster-cols TRUE
```

### Common options

- `--quant-file`: override the inferred MZmine quant CSV
- `--metadata-file`: override the inferred treated metadata TSV
- `--sample-type`: keep a specific `sample_type`, or use `all`
- `--sample-annotation`: metadata column used for the top sample annotation bar
- `--top-n`: number of annotated features to plot; use `0` for all joined features
- `--transform`: `log10` or `none`
- `--scale`: `row_zscore` or `none`
- `--cluster-rows`, `--cluster-cols`: toggle hierarchical clustering
- `--width`, `--height`: base widget size in pixels

## CANOPUS NPC Clustergrammer

Use `src/canopus_npc_clustergrammer.py` to generate a Clustergrammer-based interactive heatmap from the same CANOPUS, MZmine, and treated metadata inputs.

This version is aimed at heavier interactive browsing:

- hierarchical clustering with dendrogram navigation
- row/column search, zoom, panning, reordering, and cropping
- row categories for NPC pathway, superclass, and class
- fixed NPC colors aligned with the `biostat_toolbox` microshades mapping

### Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-helpers.txt
```

### Basic usage

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv
```

By default this command:

- infers the MZmine quant table as `results/mzmine/<batch>_quant.csv`
- infers the metadata table as `metadata/treated/<batch>_metadata.tsv`
- keeps only `sample_type == sample`
- normalizes samples first using median-sum scaling computed across all quantified features
- keeps only features where CANOPUS pathway, superclass, and class probabilities are each at least `0.85`
- ranks features by discrimination across the selected sample annotation and keeps the top 300 by default
- applies `log10(x + 1)` followed by row z-scoring
- writes `canopus_npc_clustergrammer.html` next to the CANOPUS file
- writes a companion TSV with the exact plotted data matrix and annotations

For a curated metabolomics-oriented preset, use:

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --preset metabolomics_sota
```

This preset switches the view to a more interpretable family-level heatmap:

- aggregate by NPC `superclass`
- median-sum sample normalization
- `log10(x + 1)` transform
- row-wise Pareto scaling
- discriminant ranking
- confidence filter at NPC probability `>= 0.85`
- minimum effect size and fold-change guards
- correlation distance with complete linkage

The exported table now includes discriminant evidence columns such as:

- `discriminant_score`: ANOVA F-like discrimination score
- `effect_size_eta2`: proportion of variance explained by the sample grouping
- `max_abs_log2_fc`: strongest pairwise group separation on the normalized abundance scale
- `p_value` and `fdr`: ANOVA p-value and Benjamini-Hochberg adjusted FDR

### Example

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --output-html /tmp/canopus_npc_clustergrammer.html \
  --top-n 500 \
  --cluster-cols true
```

### Common options

- `--quant-file`: override the inferred MZmine quant CSV
- `--metadata-file`: override the inferred treated metadata TSV
- `--sample-type`: keep a specific `sample_type`, or use `all`
- `--sample-annotation`: metadata column used for the main sample category
- `--preset`: `default`, `metabolomics_sota`, or `metabolomics_sota_feature`
- `--top-n`: number of annotated features to plot; use `0` for all joined features
- `--top-by`: feature ranking mode: `discriminant`, `intensity`, or `variance`
- `--aggregate-by`: aggregate rows before clustering: `none`, `class`, `superclass`, or `pathway`
- `--aggregate-method`: family aggregation summary: `sum`, `mean`, or `median`
- `--sample-normalization`: sample-wise normalization mode: `median_sum` or `none`
- `--transform`: `log10` or `none`
- `--scale`: `row_zscore`, `pareto`, or `none`
- `--dist-type`: clustering distance metric, for example `cosine`, `correlation`, or `euclidean`
- `--linkage-type`: clustering linkage, for example `average`, `complete`, or `ward`
- `--npc-prob-threshold`: minimum required CANOPUS probability at all three NPC levels; set `0` to disable
- `--min-effect-size`: minimum eta-squared effect size
- `--min-abs-log2-fc`: minimum strongest pairwise group fold-change
- `--max-pvalue`: maximum ANOVA p-value
- `--max-fdr`: maximum Benjamini-Hochberg FDR
- `--cluster-rows`, `--cluster-cols`: use clustered or alphabetical initial ordering

### Suggested discovery presets

Feature-level shortlist:

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --min-effect-size 0.15 \
  --min-abs-log2-fc 0.8 \
  --max-fdr 0.2 \
  --top-n 150
```

Superclass-level family overview:

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --aggregate-by superclass \
  --aggregate-method sum \
  --top-by discriminant
```

Metabolomics-oriented preset:

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --preset metabolomics_sota
```

Pathway-level broad overview:

```bash
python3 src/canopus_npc_clustergrammer.py \
  --canopus-file /path/to/results/sirius/canopus_structure_summary.tsv \
  --aggregate-by pathway \
  --aggregate-method sum
```

### Note

- The generated HTML loads Clustergrammer and its JavaScript dependencies from public CDNs by default, so opening the report requires network access unless you override the asset URLs.

## Spectral Modules

Use `src/analyze_spectral_modules.R` to detect groups of spectrally related features from the GNPS molecular network and test whether those modules discriminate the target factor defined in `params.yaml`.

### Basic usage

```bash
Rscript src/analyze_spectral_modules.R \
  --hash 1c63fd85afc717643ca8186ff99e6ae4
```

This command:

- reads `filtered_pairs.tsv` from the GNPS job specified by `gnps_job_id` in `params.yaml`
- builds a cosine-weighted spectral network using only features present in the selected result directory
- detects spectral modules with Louvain clustering on the network, with automatic fallback to GNPS connected components when needed
- computes a module score for each sample using the first principal component of the module feature matrix
- tests each module against the target factor with permutation p-values and reports FDR-adjusted q-values

### Common options

```bash
Rscript src/analyze_spectral_modules.R \
  --hash 1c63fd85afc717643ca8186ff99e6ae4 \
  --n-permutations 999 \
  --min-module-size 3 \
  --output-dir /tmp/spectral_modules
```

- `--hash`: hash of the result directory to read under the configured stats directory
- `--results-dir`: explicitly point to a result directory containing `DE.rds` and `foldchange_pvalues.csv`
- `--network-file`: override the default GNPS `filtered_pairs.tsv` path
- `--module-mode`: `neutral` or `discriminant`
- `--module-method`: `auto`, `louvain`, or `component`
- `--discriminant-stat`: for discriminant mode, feature-level score used to reweight spectral edges; `t` uses the raw sample matrix and `signed_log10p` uses the fold-change table
- `--same-direction-only`: for discriminant mode, optionally keep only edges connecting features with the same effect direction
- `--min-module-size`: minimum number of features required to retain a module
- `--min-cosine`: filter out network edges below a cosine threshold before module detection
- `--n-permutations`: number of label permutations used for module p-values
- `--top-n-plots`: number of top-ranked modules exported as score boxplots, GraphML networks, and PNG network previews
- `--output-dir`: directory for the module summary, member table, score table, summary plot, GraphML exports, and network previews

### Outputs

- `spectral_module_summary.csv`: one row per spectral module with module size, edge density, permutation p-value, q-value, effect size, hub feature, and dominant annotations
- `spectral_module_members.csv`: feature-to-module assignments with weighted degree, discriminant score, module loading, per-feature p-value, and annotation columns
- `spectral_module_scores.csv`: per-sample module scores used for statistical testing
- `spectral_module_scores_top.png`: boxplots of the top-ranking modules across the target groups
- `graphml/*.graphml`: one GraphML file per top-ranked module with node fold changes, discriminant scores, p-values, degree, loadings, annotations, edge cosine values, clustering weights, and exported `x/y` layout coordinates
- `network_plots/*.png`: quick network previews for the same top-ranked modules

### Notes

- `--module-mode neutral` is the baseline analysis: the network topology defines the modules, and discrimination is evaluated afterward.
- `--module-mode discriminant` is exploratory: spectral edges are reweighted by feature discrimination before clustering, so the discovered modules are explicitly biased toward connected and phenotype-associated subnetworks.
- Because discriminant mode uses phenotype information during module discovery, its downstream module p-values should be treated as exploratory rather than strict confirmatory inference.

## Optional Python helpers

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-helpers.txt
python3 src/chem_mapper.py --help
python3 src/enrich_ik.py --help
```

## Maintainer notes

Dependencies are managed with [renv](https://rstudio.github.io/renv/). The `renv.lock` file is the single source of truth — never edit it manually.

### Adding or updating a package

```bash
# CRAN package
Rscript -e "renv::install('packagename'); renv::snapshot()"

# GitHub package (pin to a specific commit for reproducibility)
Rscript -e "renv::install('user/repo@commitsha'); renv::snapshot()"

git add renv.lock
git commit -m "add packagename"
```

### Rebuilding the environment from scratch

```bash
Rscript -e "renv::restore()"
```
