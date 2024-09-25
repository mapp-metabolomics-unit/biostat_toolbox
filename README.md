# biostat_toolbox
The MAPP Biostat toolbox


## How to install the biostat_toolbox

1. Clone the current repository to your local machine

```bash
git clone https://github.com/mapp-metabolomics-unit/biostat_toolbox.git
```

2. Install the biostat_toolbox package

Here we use `renv` ton manage the libraries required for bioastat_toolbox. 

The first step once you have navigated to the cloned repository is to activate a R session.
Radian is recommended for this purpose. 

```bash
radian
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


1. Fetch the following demon repo

```bash
git clone https://github.com/mapp-metabolomics-unit/johnny-watanabe-group.git
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


