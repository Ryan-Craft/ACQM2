

one_s(x) = 2*x*exp(-x)
two_s(x) = (x/sqrt(2))*(1-(x/2))*exp(-x/2)
three_s(x) = (2*x/sqrt(27))*(1-(2*x/3)+(2*x**2/27))*exp(-x/3)

set xrange [-0.1:20]
set yrange [-0.75:0.8]
set grid

set title "l=0,α=1, 1s, 2s and 3s Hydrogen Radial Wavefunctions\n vs N=500 Laguerre Basis Approximations"
set title font ",12"

set xlabel "Radius (a_{0})" 
set ylabel "Ψ_{i}(r)"

set xlabel font ",14"
set ylabel font ",14"
set tics font ",12"
set key font ",12"


plot one_s(x) lc 1, two_s(x) lc 2, three_s(x) lc 6, 'wfout.txt' using 1:2 t "1s from program" pointtype 2 lc 8,'wfout.txt' using 1:3 t "2s from program" pointtype 6 lc 8, 'wfout.txt' using 1:4 t "3s from program" pointtype 4 lc 8


