# biostat_toolbox

Script-based metabolomics analysis utilities for the MAPP biostat workflow.

## Supported platforms

`biostat_toolbox` is currently supported on `macOS` and `Linux`.
Windows is not part of the supported installation path in this repository.

## Installation

### 1. Install system prerequisites

- `R 4.2.x`
- A working compiler toolchain for R packages
- `pandoc` available on `PATH`

### 2. Clone the repository

```bash
git clone https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
cd biostat_toolbox
```

### 3. Restore the project environment

`renv.lock` is the single source of truth for the R environment.

```bash
Rscript install.R
```

### 4. Copy the parameter templates

```bash
cp params/params_template.yaml params/params.yaml
cp params/params_user_template.yaml params/params_user.yaml
```

Edit the copied files to point to your input and output directories.

### 5. Run the preflight check

```bash
Rscript scripts/check_install.R
```

### 6. Run the toolbox

From the repository root:

```bash
Rscript src/biostat_toolbox.r
```

## Optional Python helpers

The Python scripts are optional and are not required for the main R workflow.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-helpers.txt
```

Available helper entrypoints:

```bash
python3 src/chem_mapper.py --help
python3 src/enrich_ik.py --help
```

## Useful checks

Verify that the plotting helper resolves correctly in the restored R environment:

```bash
Rscript src/plot_selected_boxplots.R --help
```

Check the current renv state:

```bash
Rscript -e 'renv::status()'
```

## Notes

- Run commands from the repository root unless a command explicitly says otherwise.
- Do not install packages from inside `src/biostat_toolbox.r`; use `Rscript install.R`.
- If the parameter files are missing, copy them again from the templates in `params/`.
