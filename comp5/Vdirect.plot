
set grid

splot 'Vdirectout.txt' u 1:2:3 with lines t "Direct Matrix"

set xrange [0:5]
set yrange [0:5]

set tics font ",12"
set title font ",12"
set xlabel font ",12"
set ylabel font ",12"
set zlabel font ",12"
set key font ",12"

set zlabel offset -2

set title "Direct Matrix 1s-3s"
set zlabel "V(k',k)" rotate by 90
set xlabel "k'"
set ylabel "k"



