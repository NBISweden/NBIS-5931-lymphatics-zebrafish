#!/usr/bin/env Rscript

# based on https://github.com/casbap/ncRNA/blob/main/docker/rpkgs.R
# used in images based on rocker/tidyverse:4.4.3

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c('TADCompare','HiCcompare'),
	version = "3.20", ask=FALSE, update=FALSE)

print("Install Bioconductor packages, done!")

