# biostat_toolbox

Script-based metabolomics analysis utilities for the MAPP biostat workflow.

## Supported platforms

`biostat_toolbox` is supported on `macOS` and `Linux`.
Windows is not part of the supported install path in this repository.

## Installation

### System prerequisites

- `R 4.2.x`
- A working compiler toolchain for R packages (Xcode CLT on macOS, `build-essential` on Linux)
- `pandoc` on `PATH`

### Steps

```bash
git clone https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
cd biostat_toolbox
Rscript install.R
cp params/params_template.yaml params/params.yaml
cp params/params_user_template.yaml params/params_user.yaml
```

Edit `params/params.yaml` and `params/params_user.yaml` for your dataset, then verify everything is in order:

```bash
Rscript scripts/check_install.R
```

### Run

```bash
Rscript src/biostat_toolbox.r
```

## Useful checks

Check the boxplot helper:

```bash
Rscript src/plot_selected_boxplots.R --help
```

Check the current renv state:

```bash
Rscript -e 'renv::status()'
```

## Optional Python helpers

The Python scripts are optional and not required for the main R workflow.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-helpers.txt
python3 src/chem_mapper.py --help
python3 src/enrich_ik.py --help
```

## Maintainer notes

### How dependencies are managed

All explicit R dependencies are declared in **`packages.yaml`** — this is the single file to edit when adding or removing a package:

```yaml
cran:
  - ggplot2
  - ...
bioc:
  - pmp
github:
  - mapp-metabolomics-unit/MAPPstructToolbox@<sha>
```

- `install.R` restores the full environment from `renv.lock` (end-user path).
- `scripts/check_install.R` reads `packages.yaml` to verify top-level packages are present.
- `scripts/rebuild_lockfile.R` reads `packages.yaml` to install everything and regenerate `renv.lock` (maintainer path).

### Updating the lockfile

After editing `packages.yaml`, regenerate `renv.lock`:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla scripts/rebuild_lockfile.R --clean
```

Use `--clean` for a full rebuild from scratch. Omit it for incremental updates.
