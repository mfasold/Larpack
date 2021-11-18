#!/usr/local/bin/gnuplot -persist
set terminal png 
set output 'PolymeraseBias.png'
set xlabel "SomLogI(total)"
set ylabel "SumLogI(transcript subset) - SumLogI(total)"
#set yrange[-0.13:0.13]
set yrange[-0.28:0.28]
#set xlabel "Probeset Average Intensity <I>"
#set ylabel "<I>(transcript subset) - <I>"
#plot "3PrimeBias.dat" u 1:2 w l title "3' End Subset", "3PrimeBias.dat" u 1:4 w l title "5' Subset", "3PrimeBias.dat" u 1:3 w l title "Transcript Middle"
# plot "3PrimeBias.dat" u 1:($8-$1) w l title "Transcript Begin", "3PrimeBias.dat" u 1:($6-$1) w l title "Transcript Middle", "3PrimeBias.dat" u 1:($3-$1) w l title "Transcript End"
plot "3PrimeBias.dat" u 1:4 w l title "Transcript End", "3PrimeBias.dat" u 1:2 w l title "Transcript Begin", "3PrimeBias.dat" u 1:3 w l title "Transcript Middle"
