#!/usr/local/bin/Rscript

# Script for installing the necessary packages
# used in DADA2manager UI.

# CRAN packages:
install.packages("optparse")

# Bioconstructor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead", version = "3.8")
BiocManager::install("dada2", version = "3.8")
BiocManager::install("phyloseq", version = "3.8")
