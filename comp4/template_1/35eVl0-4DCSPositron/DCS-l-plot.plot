
set grid

set tics font ",12"
set title font ",12"

set xlabel font ",12"
set ylabel font ",12"
set xlabel "Scattering Angle (degrees)"
set ylabel "Differential Cross Section (a_{0}^2/ Sr)"
set title "Differential Cross Section (a_{0}^2/ Sr) of Increasing l_{max}"

set key font ",12"

plot "DCSoutl0.txt" u 1:2 t "l_{max}=0", "DCSoutl1.txt" u 1:2 t "l_{max}=1","DCSoutl2.txt" u 1:2 t "l_{max}=2","DCSoutl3.txt" u 1:2 t "l_{max}=3", "DCSoutl4.txt" u 1:2 t "l_{max}=4"
