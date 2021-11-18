#!/bin/bash

# Creates eps files of all gnuplot plots ending with .gnuplot,
# Creates png files for all eps files,
# Creates graphics.pdf which saves all graphis
# (All files are saved in the subdirectory)

# Print error message
if [ $# -lt 1 ]
then
    echo "Usage: `basename $0` directory_with_arrayCorrect_output"
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
sed -e "s/HookplotName/Hookplot-Primary/ig" ${scriptDir}/HookcurvePlotTemplate.gnuplot > Hookplot-Primary.tmp1
#replace probesetStatisticsColumns with (0.5*($3+$4)):($3-$4)
sed -e "s/probesetAverages/(0.5*(\$3+\$4)):(\$3-\$4)/ig" Hookplot-Primary.tmp1 > Hookplot-Primary.tmp2
sed -e "s/smoothedHookcurve/7:8/ig" Hookplot-Primary.tmp2 > Hookplot-Primary.tmp3
sed -e "s/TheoreticCurvePlottingcommands//ig" Hookplot-Primary.tmp3 > Hookplot-Primary.tmp4

# Same for the corrected curve:
sed -e "s/HookplotName/Hookplot-Corrected/ig" ${scriptDir}/HookcurvePlotTemplate.gnuplot > Hookplot-Corrected.tmp1
sed -e "s/probesetAverages/(0.5*(\$5+\$6)):(\$5-\$6)/ig" Hookplot-Corrected.tmp1 > Hookplot-Corrected.tmp2
sed -e "s/smoothedHookcurve/9:10/ig" Hookplot-Corrected.tmp2 > Hookplot-Corrected.tmp3
if [ -f "theoreticCurve.dat" ] # Insert theoretic curve only if file exists 
then
    sed -e "s/TheoreticCurvePlottingcommands/, \"theoreticCurve.dat\" u 3:4 title \"Theoretic Curve\" w l lw 3/ig" Hookplot-Corrected.tmp3 > Hookplot-Corrected.tmp4
else
    sed -e "s/TheoreticCurvePlottingcommands//ig" Hookplot-Corrected.tmp3 > Hookplot-Corrected.tmp4
fi

${scriptDir}/replaceWithCompleteDictionary.py dataLogger.log Hookplot-Corrected.tmp4 >Hookplot-Corrected.tmp5
${scriptDir}/replaceWithCompleteDictionary.py parameters.log Hookplot-Corrected.tmp5 >Hookplot-Corrected.gnuplot
${scriptDir}/replaceWithCompleteDictionary.py dataLogger.log Hookplot-Primary.tmp4 >Hookplot-Primary.tmp5
${scriptDir}/replaceWithCompleteDictionary.py parameters.log Hookplot-Primary.tmp5 >Hookplot-Primary.gnuplot

rm Hookplot-Primary.tmp*
rm Hookplot-Corrected.tmp*

# Create images
for SCRIPTFILE in *.gnuplot
do
    gnuplot $SCRIPTFILE
done
if [ -f "3PrimeBias.dat" ] # Create transcript bias image only if file exists 
then
    echo "Creatung Tongs-plot"
    ${scriptDir}/makeZange.gnuplot
fi
mv SensitivityProfileNsPm-Corrected-single0.png SensitivityProfileNsPm-Corrected-single.png # exception for g-stacks images
mv SensitivityProfileNsMm-Corrected-single0.png SensitivityProfileNsMm-Corrected-single.png # exception for g-stacks images
mv *.png $DIRNAME  # ...and move them to graphics subdir

# Create tex file
cd $DIRNAME
${scriptDir}/replaceWithCompleteDictionary.py ../dataLogger.log ${scriptDir}/reportTemplate.tex >HookReport.tmp
${scriptDir}/replaceWithCompleteDictionary.py ../parameters.log HookReport.tmp >HookReport.tex

# remove comments, replace tab with ": ", introduce latex line-breaks (//) replace and
# mark files
cat ../parameters.log | sed '/^\#/d' | ~/projects/Larpack/scripts/replaceWithCompleteDictionary.py  ~/projects/Larpack/scripts/parameterDesciption.txt | sed -e 's|\t|: \\\\ \\> |' -e 's_> \(.*[/\\].*\)$_>  \\verb|\1|_' -e 's|$| \\\\|'  >parameters.tex

# ${scriptDir}/replaceWithCompleteDictionary.py ${scriptDir}/parameterDesciption.txt ../parameters.log HookReport.tmp >HookReport.tex
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
#rm *.tex
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