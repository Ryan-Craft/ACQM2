
set grid
set style fill solid

splot 'Exchange.txt' u 1:2:3 with lines
set xyplane 0.1

set xrange [0:5]
set yrange [0:5]

set tics font ",12"
set title font ",12"
set xlabel font ",12"
set ylabel font ",12"
set key font ",12"

set title "Exchange Matrix 1s-1s for 10 eV Projectile"
set ylabel "V(k',k)"
