#! /usr/bin/Rscript
#
# Converts celfiles into V4 format
#
# USAGE:
#  converCelfile.R ../GSE210/*.CEL
args <- commandArgs(TRUE)

library(affxparser)

for (celfile in Sys.glob(args)) {
  if (is.na(file.info(basename(celfile))$size)) {
    print(celfile)
    convertCel(celfile, basename(celfile), verbose=F)
    ## convertCel(celfile, paste(substr(basename(celfile), 1, nchar(basename(celfile))-4), ".V4.CEL", sep=""))
  } else {
    cat("File ", basename(celfile), " exists: skipping. \n")
  }
}
