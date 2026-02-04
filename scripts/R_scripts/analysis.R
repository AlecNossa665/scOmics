args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

library(Seurat)
library(dplyr)

metadata <- read.delim("/gpfs/data/courses/single_cell_omics/hrr_minimal_metadata.tsv")

