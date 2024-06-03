
set grid
set logscale y

set title font ",12"
set tics font ",11"

set xlabel font ",12"
set ylabel font ",12"

set title "Partial Integrated Cross Section (a_0^2) for Given l \nvs Incident Positron Energy (eV)"

set format y "10^{%T}"
set lmargin 10

set xlabel "Incident Energy eV"
set ylabel "Integrated Cross Section (a_0^2)" offset -1.5

set xrange [0:50]
set yrange [0.00001:100]

plot "collated_ICS.txt" u 1:2 w line t "l=0", "collated_ICS.txt" u 1:3 w line t "l=1", "collated_ICS.txt" u 1:4 w line t "l=2", "collated_ICS.txt" u 1:5 w line t "l=3"

