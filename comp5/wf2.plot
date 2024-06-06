

two_p(x) = (x**2/sqrt(24))*exp(-x/2)
three_p(x) = (8*x**2/(27*sqrt(6)))*(1-(x/6))*exp(-x/3)
four_p(x) = (x/(64*sqrt(15))) * x * (0.25*x**2 - 5*x +20)*exp(-x/4)


set xrange [-0.1:35]
set yrange [-0.75:0.8]
set grid

set title "l=1,α=1, 2p, 3p and 4p Hydrogen Radial Wavefunctions\n vs N=500 Laguerre Basis Approximations"
set title font ",12"

set xlabel "Radius (a_{0})"
set ylabel "Ψ_{i}(r)"

set xlabel font ",14"
set ylabel font ",14"
set tics font ",12"
set key font ",12"


plot two_p(x), three_p(x), four_p(x) , 'wfout.txt' using 1:2 t "2p from program", 'wfout.txt' using 1:3 t "3p from program", 'wfout.txt' using 1:4 t "4p from program"


