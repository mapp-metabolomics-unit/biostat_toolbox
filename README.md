# biostat_toolbox
The MAPP Biostat toolbox


## How to install the biostat_toolbox

Note this has been tested on UNIX platforms (Ubuntu and MacOS non ARM).
Not guaranteed to work in other setups.


### Clone the current repository to your local machine

```bash
git clone https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
```
### Install the cappropriate R version (4.2.2)

For this we recommend using [rig](https://github.com/r-lib/rig) 

### Install the biostat_toolbox package

Here we use [renv](https://github.com/rstudio/renv) to manage the libraries required for bioastat_toolbox. 

```bash
install.packages("renv")
```

The first step once you have navigated to the cloned repository is to activate a R session.
Radian is recommended for this purpose. 

```bash
radian
```

Alternatively R works of course

```bash
R
```

Then, you can activate the `renv` environment by running the following command:

```R
renv::activate()
```

Once the `renv` environment is activated, you can install the `biostat_toolbox` package by running the following command:

```R
renv::restore()
```

This should install all the required libraries for the `biostat_toolbox` package.


If you are having trouble with the installation, please drop an issue [here](https://github.com/mapp-metabolomics-unit/biostat_toolbox/issues)


## Example usage


1. Fetch the following demo repo

```bash
git clone https://github.com/mapp-metabolomics-unit/johnny-watanabe-group.git
```

or, if you use ssh authentification.

```bash
git clone git@github.com:mapp-metabolomics-unit/johnny-watanabe-group.git
```

2. Run the MAPP copier template `update` function to update paths.

```bash
cd johnny-watanabe-group
copier update --trust
```

You should now be able to use the updated https://github.com/mapp-metabolomics-unit/johnny-watanabe-group/blob/46c665c607fc3137f5f2de65288bcda11ab13a93/docs/mapp_project_00051/mapp_batch_00114/CONTRIB.md#run-biostat_toolbox params and run biostat_toolbox functions.

```bash
Rscript src/biostat_toolbox.R
```




