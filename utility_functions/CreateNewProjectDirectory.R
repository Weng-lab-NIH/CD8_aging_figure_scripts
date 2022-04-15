#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
# DESCRIPTION:
# This Script creates a new directory structure when beginning a new project.

library(tidyverse)

path <- as.character(args[1])

folder <- function (name) {
  dir.create(paste(path, name, sep = '/'))
}

# CORE FOLDER STRUCTURE (DO NOT CHANGE THESE!!!)

# Data Files

folder('data')
folder('data/raw')
folder('data/tidy')
folder('data/interim')
folder('data/external')

# Documentation

folder('docs')
folder('docs/notebooks')
folder('docs/references')
folder('docs/paper')
folder('docs/presentations')

# Models

folder('models')

# Code

folder('code')
folder('code/raw_code')
folder('code/final_code')

# Figures

folder('figures')
folder('figures/exploratory_figures')
folder('figures/tidy_figures')


# OPTIONAL FOLDERS:

folder('metadata')