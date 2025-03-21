#!/usr/bin/env Rscript

#' Package management module for biostat toolbox
#' This module handles all package installations and loading

#' Helper function to load or install packages
#' @param p Package name to install/load
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {
    install.packages(p, dep = TRUE, Ncpus = 40)
  }
  require(p, character.only = TRUE)
}

#' Initialize all required packages for the biostat toolbox
initialize_packages <- function() {
  # Set default CRAN repo
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.us.r-project.org"
  options(repos = r)
  rm(r)
  
  # Install BiocManager if not present
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Required packages organized alphabetically
  required_packages <- c(
    "crosstalk",
    "digest",
    "dplyr",
    "DT",
    "ggh4x",
    "ggrepel",
    "htmltools",
    "iheatmapr",
    "janitor",
    "microshades",
    "plotly",
    "pls",
    "pmp",
    "readr",
    "rfPermute",
    "tidyverse",
    "vegan",
    "webchem",
    "wesanderson",
    "WikidataQueryServiceR",
    "yaml"
  )
  
  # Load all required packages
  sapply(required_packages, usePackage)
  
  # Load MAPPstructToolbox
  library(MAPPstructToolbox)
} 