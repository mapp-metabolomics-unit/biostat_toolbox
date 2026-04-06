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
