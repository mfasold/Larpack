## ---------------------------------------------------------
## Method for normalising GeneChip gene expression data using hook curves
## (C) Joerg Hackermueller, Mario Fasold and Jan Bruecker
## ---------------------------------------------------------

## Usage:
##   library(larpack)
##   (1) Using Affybatch object:
##   library(affy)
##   data = ReadAffy()
##   l = runLarpack(data[,1], "HG_probeTab", list("probeset-limit"="8000")) 
## 
##   (2) Direct Import
##   l = runLarpack("ChipA.cel", "HG_probeTab", list("probeset-limit"="8000"))

runLarpack <- function (celFileOrAffybatch,
                         seqFile,
                         optionsList=list(),
                         .checkArgs=TRUE)
{
  useAffybatch <- FALSE 
  
  if (.checkArgs) {
    ## Check if affybatchObject provided
    if (is.object(celFileOrAffybatch)) {
      useAffybatch <- TRUE
      ## @todo Insert tests
    }
    else {  ## No affybatch object -> must read celfile      
      ## Check whether celFileOrAffybatch is OK
      if(length(celFileOrAffybatch) != 1) {
        stop("Argument 'celFile' should be a single file: ", 
             paste(celFileOrAffybatch, collapse=", "));
      }
      ## Expand '~' pathnames to full pathnames.
      celFileOrAffybatch <- file.path(dirname(celFileOrAffybatch), basename(celFileOrAffybatch));
      if (!file.exists(celFileOrAffybatch)) {
        stop("Cannot read signal intensities. Files not found: ", celFileOrAffybatch);
      }

    }
    ## Check whether seqFile is OK
    if(length(seqFile) != 1) {
      stop("Argument 'seqFile' should be a single file: ", 
           paste(seqFile, collapse=", "));
    }
    ## Expand '~' pathnames to full pathnames.
    seqFile <- file.path(dirname(seqFile), basename(seqFile));
    if (!file.exists(seqFile)) {
      stop("Cannot read sequence db file. File not found: ", seqFile);
    }
    
  } # if (.checkArgs)

  if (useAffybatch) {
    chipCount = length(celFileOrAffybatch)

    # Create the return variable (list) and get the matrix of intensities, which are updated
    larDat = list(affybatch = celFileOrAffybatch)
    probeExpressionMeasures = exprs(celFileOrAffybatch)

    for(chipIndex in seq(1, chipCount)) {
      sampleName = sampleNames(celFileOrAffybatch)[chipIndex]
      cat("Processing chip ", sampleName, " (", chipIndex ,"/", chipCount,").\n")
      
      singleDat <- .Call("larpackOligo",
                         "/useAffybatch",
                         seqFile,
                         optionsList,
                         as.vector(exprs(celFileOrAffybatch[, chipIndex])),
                         nrow(celFileOrAffybatch), ncol(celFileOrAffybatch), 
                         PACKAGE="larpack");

      ## Update intensities with expression measures
      probeExpressionMeasures[,chipIndex] = singleDat$ProbeExpressionMeasures

      ## ## Write larpack statistics into list element `[chipIndex]`
      ## larDat[[as.character(chipIndex)]] = singleDat

      ## Create data frame for each of the lists returned by larpack
      excludeElements = c("ProbesetId", "ProbeExpressionMeasures")
      elements = setdiff(names(singleDat),  excludeElements)
      for (el in elements) {
        if (chipIndex == 1) { # Create data frame for list elements in first iteration                
          larDat[[el]] = data.frame(singleDat[[el]])
          colnames(larDat[[el]]) = sampleName
          if (length(singleDat[[el]]) == length(singleDat[["ProbesetId"]])) { # assign probeset names if possible
            rownames(larDat[[el]]) = singleDat[["ProbesetId"]]
          }
        } else { # append column for that list element
          larDat[[el]][sampleName] = singleDat[[el]]
        }
      }
    }

    ## Update all expression measures
    exprs(larDat$affybatch) = probeExpressionMeasures
  }
  else {
    larDat <- .Call("larpackOligo",
                    celFileOrAffybatch,
                    seqFile,
                    optionsList,
                    NA,0,0,
                    PACKAGE="larpack");
  }
  return (larDat);
  
} # runLarpack  


## Computes the threshold separating non-specific ans specific binding
## using the "change point of the hook curve". 
##
computeNsThreshold <- function (x,y, method = 2, digitizeIntervals = 100)
{
  larDat <- .Call("larpackComputeNsThreshold",
                  x,y, as.integer(method),
                  as.integer(digitizeIntervals),
                  PACKAGE="larpack");
  return (larDat);
  
} # computeNsThreshold


computeSequenceProfile <- function (sequences, probesets, intensities, rank = 2)
{
  larDat <- .Call("larpackComputeSequenceProfile",
                  sequences, probesets, intensities, as.integer(rank), 
                  PACKAGE="larpack");
  return (larDat);
  
} # computeSequenceProfile

