# biostat_toolbox

Script-based metabolomics analysis utilities for the MAPP biostat workflow.

## Supported platforms

`biostat_toolbox` is supported on `macOS` and `Linux`.
Windows is not part of the supported install path in this repository.

## Simple installation

Install these system prerequisites first:

- `R 4.2.x`
- a working compiler toolchain for R packages
- `pandoc` on `PATH`

Then from a fresh clone:

```bash
git clone https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
cd biostat_toolbox
Rscript install.R
cp params/params_template.yaml params/params.yaml
cp params/params_user_template.yaml params/params_user.yaml
Rscript scripts/check_install.R
```

Edit `params/params.yaml` and `params/params_user.yaml` for your machine, then run:

```bash
Rscript src/biostat_toolbox.r
```

## Useful checks

Check that the helper script resolves correctly:

```bash
Rscript src/plot_selected_boxplots.R --help
```

Check the current `renv` state:

```bash
Rscript -e 'renv::status()'
```

## Optional Python helpers

The Python scripts are optional and are not required for the main R workflow.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-helpers.txt
python3 src/chem_mapper.py --help
python3 src/enrich_ik.py --help
```

## Maintainer notes

- `renv.lock` is the single source of truth for the R environment.
- This project disables `pak` during `renv` operations.
- If you intentionally change R dependencies, rebuild the lockfile with:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla scripts/rebuild_lockfile.R --clean
```
