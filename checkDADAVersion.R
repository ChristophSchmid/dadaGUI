#!/usr/local/bin/Rscript

#Short script for checking the version of the installed dada2 package

if("dada2" %in% rownames(installed.packages())) {
  library(dada2)
  write(unname(getNamespaceVersion("dada2")), stdout())
} else {
  write("none installed", stdout())
}
