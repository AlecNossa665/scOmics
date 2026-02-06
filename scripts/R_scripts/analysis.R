args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

library(Seurat)
library(dplyr)

metadata <- read.delim("/gpfs/data/courses/single_cell_omics/hrr_minimal_metadata.tsv")
# TODO create a data frame with like patient id, sample, tissue type

# we can use an apply function, which takes a list of vectors and applies a function to each item in the vector.
# treat the metadata as a matrix and use indexing notation
# barcodes are not unique across libraries, we have to rename each cell to something that is unqiue. Use the RenameCell object to construct. 
# new_names = pastr(row_names(serut_metadata_use)[i]
#                                           colnames(seur_obj), sep = '_'

# when you have more than 2 things to collapse use the R reduce function and then use the seurat function merge

#work flow
# load object, which is individual smaples
# merge meta data
# make unique, because barcodes duplicate between samples 
# merge all the objects with the unique bar codes

#seurat


