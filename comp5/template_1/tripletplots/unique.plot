set grid

set title font ",12"
set ylabel font ",12"
set xlabel font ",12"
set tics font ",12"
set key font ",11"

set ylabel "Triplet 1s-1s Matrix Elements"
set xlabel "Off shell linear momentum k (a.u.)"

set pointsize 1.5

set yrange [-1.5:0.5]

set title "Half-Offshell V and K matrix elements for Kgrid 1 and 2"

plot "Koff-grid1th1.txt" u 1:2 w lp lc 1 pointtype 4 t "kgrid 1 K(k), θ=1", "Vmat-grid1th1.txt" u 1:2 w lp lc 2 pointtype 4 t "kgrid 1 V(k) θ=1", "Koff-grid2th2.txt" u 1:2 w lp lc 1 pointtype 8 t "kgrid 2 K(k) θ=2", "Vmat-grid2th2.txt" u 1:2 w lp lc 2 pointtype 8 t "kgrid 2 V(k) θ=2"


