# The Leipzig Array Package - Installation Guide

Author: Mario Fasold

## Introduction
Microarray signal intensities are affected by various
interferences. The Larpack program package aims at correcting these
intensities to obtain better measures for target abundance. The method
is based on a physico-chemical model that incorporates
cross-hybridization, sequence effects and saturation. The Larpack
currently supports Affymetrix (c) 3' expression (genechip) and tiling
microarray data that contains perfect-match and mismatch probe pairs.

## Downloading and Using the Windows Binaries
The program does not need to be installed. Instead, the executable is
stored in an arbitrary directory:
  1. Download and install the [Visual C++ Redistributable](http://www.microsoft.com/downloads/details.aspx%3FFamilyId%3D32BC1BEE-A3F9-4C13-9C99-220B62A191EE&displaylang%3Den)
  2. Download the Larpack.zip from the [Hook-Website](https://web.archive.org/web/20161011085215/http://www.izbi.uni-leipzig.de/downloads_links/programs/hook.php) and extract it, for
     example in `C:\`. You should now have a directory `C:\Larpack`
     containing the files.
  3. Optionally add the above directory to your PATH environment
     variable, so it can be started from any directory with
     `arrayCorrect.exe` (as described [here](http://vlaurie.com/computers2/Articles/environment.htm)).

Invoke the program as follows:
  1. Create a directory to contain the results, e.g. `C:\results\chipA`.
  2. Open the Windows-Command-line with Start -> Run -> `cmd` -> OK
  3. Go to your results directory, e.g. `cd C:\results\chipA`
  4. Assuming your Celfile is `C:\data\chipA.CEL`, the
     probesequence-file for your chip is `C:\data\HG-U133A.probe_tab`
     and Larpack was extracted to `C:\Larpack`, use the following
     command to start Larpack:
     ```
     C:\Larpack\arrayCorrect.exe -i C:\data\chipA.CEL -c C:\data\HG-U133A.probe_tab
     ```
  5. (The script forAllCelfiles.bat can be used to process many
     chips of a series.)
     
## Building from Source (Linux)
To run the program under Linux you'll have to build it from source.
  1. Download the source tarball and unpack it, for example using the
     command `tar -xvzf Larpack.tar.gz`
  2. Install the following packages including their dependencies:
      - gsl-devel
      - boost-devel
       - (R-devel) (only needed for R integration)
  3. Compile the program by running `make` in the directory you
     exctracted the source (requires g++ and make). This creates the
     `arrayCorrect` program for expression and tiling arrays. Invoke
     `make pmOnlyCorrect` for the PM-only chips such as exon or Gene
     ST arrays.

## Installing Larpack as a R package
The R package is essentially a wrapper for the `arrayCorrect` program
and returns the results in appropriate data structures, such as an
`AffyBatch` object (Bioconductor package Affy). You can build and
install the R package under Linux as desribed in the file
`install.txt` in the `R` subdirectory of the source. There are no
Windows binaries yet.

## Further documentation
There is a manpage (linked at the [Hook-Website](https://web.archive.org/web/20161011085215/http://www.izbi.uni-leipzig.de/downloads_links/programs/hook.php) that describes the
program usage and available program options. Entering "make doc" in
the source directory generates a developer documentation using doxygen
(requires doxygen package).
