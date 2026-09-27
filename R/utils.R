# R/utils.R - Shared utility functions

load_bio_packages <- function() {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(BiocManager)
  })
  message('Core bioinformatics libraries loaded.')
}
