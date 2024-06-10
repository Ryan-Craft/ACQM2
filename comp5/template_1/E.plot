
set grid
set style fill solid

splot 'Exchange.txt' u 1:2:3 with lines t "Exchange Matrix 1 eV"
set xyplane 0.1

set xrange [0:5]
set yrange [0:5]

set tics font ",12"
set title font ",12"
set xlabel font ",12"
set ylabel font ",12"
set zlabel font ",12" rotate by 90
set key font ",12"

set zlabel offset -2

set title "Exchange Matrix 1s-3s for 1 eV Projectile"
set zlabel "V(k',k)"
set xlabel "k'"
set ylabel "k"
