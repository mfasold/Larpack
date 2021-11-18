set term png large size 800,600
set output "HookplotName.png" 
# set term epslatex color "Helvetica" 8
# set output "HookplotName.eps"
set title 'HookplotName IP: HookplotName-intersectionPointXvalue'
set grid lt 0 lw 1
show grid
set xrange[0.5:4.5]
set yrange[-0.2:1.2]
set xlabel '0.5*(log(PM) + log(MM))'
set ylabel 'log(PM) - log(MM)'
set key box
plot "probesetStatistics" using probesetAverages w dots title 'raw plot', "probesetStatistics" using smoothedHookcurve w points pt 3 title 'averaged curve' HookplotName-plottingCommands TheoreticCurvePlottingcommands