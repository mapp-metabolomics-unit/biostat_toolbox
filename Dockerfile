# This Dockerfile is used to build an image of the biostat_toolbox.R which can be used across multiple platform.
# It will allow the installation of the following packages in R:

# install.packages("cowplot", Ncpus = 40)
# install.packages("crosstalk", Ncpus = 40)
# install.packages("data.table", Ncpus = 40)
# install.packages("datasets", Ncpus = 40)
# install.packages("dbscan", Ncpus = 40)
# install.packages("digest", Ncpus = 40)
# install.packages("dplyr", Ncpus = 40)
# install.packages("DT", Ncpus = 40)
# install.packages("emmeans", Ncpus = 40)
# install.packages("EnhancedVolcano", Ncpus = 40)
# install.packages("fpc", Ncpus = 40)
# install.packages("ggdendro", Ncpus = 40)
# install.packages("ggh4x", Ncpus = 40)
# install.packages("ggplot2", Ncpus = 40)
# install.packages("ggraph", Ncpus = 40)
# install.packages("ggrepel", Ncpus = 40)
# install.packages("ggtree", Ncpus = 40)
# install.packages("graphlayouts", Ncpus = 40)
# install.packages("gridExtra", Ncpus = 40)
# install.packages("gt", Ncpus = 40)
# install.packages("heatmaply", Ncpus = 40)
# install.packages("here", Ncpus = 40)
# install.packages("htmltools", Ncpus = 40)
# install.packages("igraph", Ncpus = 40)
# install.packages("iheatmapr", Ncpus = 40)
# install.packages("janitor", Ncpus = 40)
# install.packages("jsonlite", Ncpus = 40)
# install.packages("magick", Ncpus = 40)
# install.packages("manhattanly", Ncpus = 40)
# install.packages("microshades", Ncpus = 40) ### remotes::install_github("KarstensLab/microshades", dependencies=TRUE)
# install.packages("modEvA", Ncpus = 40)
# install.packages("openxlsx", Ncpus = 40)
# install.packages("plotly", Ncpus = 40)
# install.packages("pls", Ncpus = 40)
# install.packages("pmp", Ncpus = 40)
# install.packages("purrr", Ncpus = 40)
# install.packages("randomcoloR", Ncpus = 40)
# install.packages("randomForest", Ncpus = 40)
# install.packages("rcdk", Ncpus = 40)
# install.packages("RColorBrewer", Ncpus = 40)
# install.packages("readr", Ncpus = 40)
# install.packages("reshape2", Ncpus = 40)
# install.packages("reticulate", Ncpus = 40)
# install.packages("rfPermute", Ncpus = 40)
# install.packages("rgl", Ncpus = 40)
# install.packages("rockchalk", Ncpus = 40)
# install.packages("ropls", Ncpus = 40)
# install.packages("this.path", Ncpus = 40)
# install.packages("tidyr", Ncpus = 40)
# install.packages("tidyverse", Ncpus = 40)
# install.packages("tools", Ncpus = 40)
# install.packages("treemap", Ncpus = 40)
# install.packages("vegan", Ncpus = 40)
# install.packages("viridis", Ncpus = 40)
# install.packages("vroom", Ncpus = 40)
# install.packages("webchem", Ncpus = 40)
# install.packages("wesanderson", Ncpus = 40)
# install.packages("WikidataQueryServiceR", Ncpus = 40)
# install.packages("yaml", Ncpus = 40)


# Use an official R 4.0 runtime as a parent image

FROM rocker/r-ver:4.0.3

# Set the working directory in the container

WORKDIR /biostat_toolbox

# Copy the current directory contents into the container at /biostat_toolbox

COPY . /biostat_toolbox

# Install any needed packages specified in the biostat_toolbox.R

RUN Rscript biostat_toolbox.R

# Make port 80 available to the world outside this container

EXPOSE 80

# Define environment variable


