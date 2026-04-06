# biostat_toolbox

Script-based metabolomics analysis utilities for the MAPP biostat workflow.

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
