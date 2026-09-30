# V2 environment

V2 uses the pinned R 4.6.1 and its own `renv.lock`. The historical root environment still pins `MAPPstructToolbox` for legacy runs; V2 uses CRAN `pls` for supervised modeling and `plotly` for interactive views, and does not load either the fork or Bioconductor `structToolbox`.

From `v2/`, install the pinned packages once, then use the same R version for CLI and viewer:

```bash
rig run -r 4.6.1 -f bootstrap.R
rig run -r 4.6.1 -f run.R -- validate --dataset ../configs/mapp_batch_00275.dataset.yaml --recipe ../configs/mapp_batch_00275.recipe.yaml
rig run -r 4.6.1 -f run.R -- plan --dataset ../configs/mapp_batch_00275.dataset.yaml --recipe ../configs/mapp_batch_00275.recipe.yaml
rig run -r 4.6.1 -f run.R -- run --dataset ../configs/mapp_batch_00275.dataset.yaml --recipe ../configs/mapp_batch_00275.recipe.yaml
```

The example `design.contrasts: all` runs 15 unordered group pairs; omitting `contrasts` also defaults to all when differential analysis is enabled. Explicit contrast entries still limit the run. Select the new run hash after updating a recipe.

To browse the completed 00275 run from `v2/`, keep this command running and open the **Listening on http://127.0.0.1:PORT** URL printed in its terminal:

```bash
MAPP_STATS_ROOT=../../didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275 \
  rig run -r 4.6.1 -e 'shiny::runApp("../app", host="127.0.0.1", launch.browser=FALSE)'
```

The root must be the batch directory or an approved parent containing MAPP project/batch folders, not `results/stats_v2` or a hash. The newest completed run is selected on startup or **Refresh catalogue**; check the hash printed by the CLI against **Completed batch / run**. **All groups** is active by default; uncheck it to select visible groups, or double-click a plot legend entry to isolate its group. PCA, PCoA and validated PLS-DA scores support independent 2D X/Y axes, plus a rotatable 3D view with a Z axis when three saved components exist. **PCA sample scores are samples**; click a point in the separate PCA *feature loadings* plot, a volcano feature or a feature ID in a table, or choose a feature in the left sidebar to open the collapsible right-side **Feature details** box plot and description. The panel narrows the central plots while open rather than obscuring them; the boxes display group medians, quartiles and whiskers, with individual sample dots. When `annotations.horizontal` was supplied, **Saved InChIKey2D consensus** filters recorded 0–3 agreeing-source counts. All source-specific `_SMILES` columns can be depicted through the public Natural Products API; exact duplicate SMILES share one card with every source label. Met-annot-enhancer structures show their reported taxon and clickable structure/taxon Wikidata links when saved. The browser sends SMILES to the external depiction service. Older runs without the horizontal table cannot offer the consensus filter until a new immutable run is published. Use the Plotly toolbar for a PNG of the current view, or the verified saved PNG/PDF buttons for immutable original 2D PCA/PCoA/volcano plots. See the [viewer walkthrough](../docs/STATS_PIPELINE_V2.md#read-only-explorer) for each tab and download steps.

**Every quantified feature remains available at the raw stage**, even if blank or QC filtering excluded it from downstream analysis. The drawer automatically labels and shows saved raw peak height/area when the selected stage lacks the feature; the newest run also exports input metadata for filtered features. Batch 00275 contains peak **height** but no peak area.

Shiny chooses a free port. On a remote server, [forward its printed port over SSH](../docs/STATS_PIPELINE_V2.md#access-from-another-machine) before opening the URL on your laptop.

The extra `--` forwards options through `rig`. Install R 4.6.1 separately if it is not already available; `Rscript run.R` is an alternative when that executable selects R 4.6.1. `renv::status()` currently warns that recorded packages are unused because it scans `v2/` while the analysis source is in `../src/v2/` and the explorer in `../app/`. Bootstrap restores and checks the required packages despite this detection limit. Keep the root lockfile separate. See [the pipeline guide](../docs/STATS_PIPELINE_V2.md) for input mappings, the approved-root viewer and scientific limits.
