# Larpack: Hook analysis of GeneChip microarrays

Status: Archived project, active between 2007 and 2016

## Description
Calibration of microarray measurements aims at removing systematic biases from the probe-level data to get expression estimates which linearly correlate with the transcript abundance in the studied samples. We address hybridization on microarrays as a reaction process in a complex environment and express the measured intensities as a function of the input quantities of the experiment.

The hook-method is a calibration approach which is based on a graphical summary of the actual hybridization characteristics of a particular microarray. The hook method provides a set of chip summary characteristics which evaluate the performance of a given hybridization. The theory and algorithm of the method is described in [this paper](http://www.almob.org/content/3/1/12).

The C++ implementation "Larpack" includes the most recent procedures producing most robust results. 

## G-stack correction

The correction of CEL-files for the effect of G-stacks (described in the [BMC Bioinformatics publication](http://www.biomedcentral.com/1471-2105/11/207)) is implemented in Larpack. Please use arrayCorrect with the command-line option "--update-celfiles".

## Releases: Larpack (C++/R implementation)

This repository contains Revision 493 of the original SVN repository (from Jun 2010). The source code contains contains arrayCorrect, pmOnlyCorrect, correctCelfile and the R Package. The windows binary arrayCorrect is available upon request, or can be build using the instructions in the repository (needs the Visual C++ Redistributable).

- Short Installation Guide ([HTML](INSTALLATION.md))
- Manpage describing program options ([HTML](https://web.archive.org/web/20161011085215/http://www.izbi.uni-leipzig.de/downloads_links/downloads/arrayCorrect.html)

## Publications

Please find further publications regarding Larpack analysis [here](https://web.archive.org/web/20161011085215/http://www.izbi.uni-leipzig.de/downloads_links/programs/hook.php)
