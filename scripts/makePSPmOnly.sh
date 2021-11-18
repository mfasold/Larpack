#!/bin/bash

# Creates eps files of all gnuplot plots ending with .gnuplot,
# Creates png files for all eps files,
# Creates HookReport.pdf which saves all graphis
# (All files are saved in the subdirectory)

# Print error message
if [ $# -lt 1 ]
then
    echo "Usage: `basename $0` directory_with_pmOnlyCorrect_output"
    exit $E_BADARGS
fi

# Get absolute path of script
initDir=`pwd` # Save current dir
cd `dirname $0` # Go to script dir 
scriptDir=`pwd` # Save script dir
cd ${initDir} # Go back 

cd $1
DIRNAME="graphics"
mkdir $DIRNAME


# Create a gnuplot script that plots the primary hookcurve. We use a template script and replace some strings:
#replace HookplotName with Hookplot-Primary:
sed -e "s/HookplotName/Hookplot-Primary/ig" ${scriptDir}/HookcurvePlotTemplatePmOnly.gnuplot > Hookplot-Primary.tmp1
sed -e "s/smoothedHookcurve/5:6/ig" Hookplot-Primary.tmp1 > Hookplot-Primary.tmp2

# Same for the corrected curve:
sed -e "s/HookplotName/Hookplot-Corrected/ig" ${scriptDir}/HookcurvePlotTemplatePmOnly.gnuplot > Hookplot-Corrected.tmp1
sed -e "s/smoothedHookcurve/7:8/ig" Hookplot-Corrected.tmp1 > Hookplot-Corrected.tmp2

${scriptDir}/replaceWithCompleteDictionary.py dataLogger.log Hookplot-Corrected.tmp2 >Hookplot-Corrected.tmp3
${scriptDir}/replaceWithCompleteDictionary.py parameters.log Hookplot-Corrected.tmp3 >Hookplot-Corrected.gnuplot
${scriptDir}/replaceWithCompleteDictionary.py dataLogger.log Hookplot-Primary.tmp2 >Hookplot-Primary.tmp3
${scriptDir}/replaceWithCompleteDictionary.py parameters.log Hookplot-Primary.tmp3 >Hookplot-Primary.gnuplot

rm Hookplot-Primary.tmp*
rm Hookplot-Corrected.tmp*

# Create images
for SCRIPTFILE in *.gnuplot
do
    gnuplot $SCRIPTFILE
done
if [ -f "3PrimeBias.dat" ] # Create transcript bias image only if file exists 
then
    ${scriptDir}/makeZange.gnuplot
fi
mv *.png $DIRNAME  # ...and move them to graphics subdir

# Create tex file
cd $DIRNAME
${scriptDir}/replaceWithCompleteDictionary.py ../dataLogger.log ${scriptDir}/reportTemplate.tex >HookReport.tmp
${scriptDir}/replaceWithCompleteDictionary.py ../parameters.log HookReport.tmp >HookReport.tex

# remove comments, replace tab with ": ", introduce latex line-breaks (//) replace and
# mark files
cat ../parameters.log | sed '/^\#/d' | ~/projects/Larpack/scripts/replaceWithCompleteDictionary.py  ~/projects/Larpack/scripts/parameterDesciption.txt | sed -e 's|\t|: \\\\ \\> |' -e 's_> \(.*[/\\].*\)$_>  \\verb|\1|_' -e 's|$| \\\\|'  >parameters.tex

rm HookReport.tmp

# Create html-file
${scriptDir}/replaceWithCompleteDictionary.py ../dataLogger.log ${scriptDir}/reportTemplate.html >HookReport.tmp
${scriptDir}/replaceWithCompleteDictionary.py ../parameters.log HookReport.tmp >HookReport.html
rm HookReport.tmp
cp ${scriptDir}/styles.css styles.css


#FILES=""
#This create png for all plots:
# for GRAPHIC in `find *.eps`
# do
# #   FILES=$FILES" "$GRAPHIC
# #   psmerge -ographics.ps $FILES
# #   ps2pdf graphics.ps graphics.pdf
# 	echo "$GRAPHIC to png"
# 	gs -dNOPAUSE -sDEVICE=png256 -sOutputFile=$GRAPHIC".png" -q -dBATCH $GRAPHIC
# done

pdflatex HookReport.tex
#latex HookReport.tex
#dvipdf HookReport.dvi
# latex2html HookReport.tex
# rm *.tex
#rm *.eps
rm *.log
rm *.aux

cd ..

# Create alternative-format profiles
mkdir profiles
for f in Sensitivity*.dat; 
do 
    ${scriptDir}/convertProfileToTableForm.py $f > profiles/$f.table
done

cd ${initDir}