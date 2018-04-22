# Run this R script with the following arguments:
# Rscript data_file link_type k 

# Clear the workspace
rm(list = ls())

# Get command line arguments
args = commandArgs(trailingOnly = TRUE)

# file = data file ((excluding begin and end))tsv) format
file = read.table(file = args[1], sep = '\t', header = FALSE)

# extract measurements of each gene
dat = as.matrix(file[,3:ncol(file)])
row.names(dat) = file[,1]

# get link type 
# S = single link
# C = complete link
# A = average link
linkType = args[2]

# get number of specified clusters that you want to have
k = as.numeric(args[3])

#####
#debug
#file = read.table(file = 'tiny-yeast.tsv', sep = '\t', header = FALSE)
#dat = as.matrix(file[,3:ncol(file)])
#row.names(dat) = file[,1]
#linkType = "S"
#k = 2

#####

source("functions.R")
result = hierarchicalClutersting(dat, file, linkType, k)
printClusterInfo(result$ListOfClusters, file)


