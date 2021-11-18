#!/bin/bash
# ---------------------------------------------------------------------------
# Enter the directory of your executable here
# ---------------------------------------------------------------------------
larpackDir="/homes/brauerei/mario/projects/Larpack"
RUN_OLIGO=$larpackDir/bin/arrayCorrect64
makePs=$larpackDir/scripts/makePS.sh

# Check if main prog exists
if [ ! -f $RUN_OLIGO ]
then
    echo "Executable does not exist. Please edit your directory"
    exit $E_BADARGS
fi

# Runs OLIGO pogram for all .cel files in a certain directory
# and stores each result in a own directory
# Print error message
if [ $# -lt 2 ]
then
    echo "Usage: `basename $0` DIRECTORY_WITH_CELFILES TABULAR_FILE ADDITIONAL_PARAMTERS"
    exit $E_BADARGS
fi

celDir=$1
tabFile=$2
shift 2 # shifts the input parameters two places left

for CELFILE in `find $celDir -maxdepth 1 -name "*.cel" -o -name "*.CEL"`
do
    # Exit when problems reading a CELFILE
    if [ ! -r $CELFILE ]
    then
        echo "Error reading file $CELFILE. Exiting."
        exit
    fi

    # Create new directory for the results of current celfile
    # DIRNAME=${file%.*}_result
    DIRNAME=${CELFILE##/*/}
    DIRNAME=${DIRNAME%.*}
    mkdir $DIRNAME
    cd $DIRNAME

    # Skip current directory expression measure table already exists
    # -> useful if computation aborted
    if [ -r ExpressionMeasure-Pm.dat ]
    then
        echo $CELFILE" has been processed already -> skipping."
        cd ..
        continue
    fi

     # Run oligo
    echo "Running $RUN_OLIGO -i $CELFILE -c $tabFile $*"
    echo " -> running in directory $DIRNAME"

    $RUN_OLIGO -i $CELFILE -c $tabFile $* >out.log

    # Create report
    $makePs .

    cd ..
done

# Create summary table containing expression measures and Hook statistics
$larpackDir/scripts/multiPaste.py 1 2 `find . -name "ExpressionMeasure-PmExp10.dat" -print |xargs ls` >allExpresssionMeasures.dat

sed -i '/^$/d' `find . -name "dataLogger.log" -print |xargs ls` # remove empty lines
$larpackDir/scripts/multiPaste.py 1 2 `find . -name "dataLogger.log" -print |xargs ls` >allDataLoggers.dat
sed -i -e '/plottingCommands/d' -e "/startingTime/d" allDataLoggers.dat # remove lines prohibiting automatic reading in R

python $larpackDir/scripts/summaryHTML.py . >summary.html