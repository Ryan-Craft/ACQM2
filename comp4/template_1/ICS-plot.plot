
set grid xtics ytics mytics
set logscale y

set title font ",12"
set tics font ",12"

set xlabel font ",12"
set ylabel font ",12"

set title "Integrated Cross Section (a_0^2) vs Incident Positron Energy (eV)\n for Different Sums of l_{max}"

set xlabel "Incident Energy eV"
set ylabel "Integrated Cross Section (a_0^2)"

set xrange [0:50]

plot "collated_ICS.txt" u 1:6 w line t "Sum to l_{max}=0", "collated_ICS.txt" u 1:7 w line t "Sum to l_{max}=1", "collated_ICS.txt" u 1:8 w line t "Sum to l_{max}=2", "collated_ICS.txt" u 1:9 w line t "Sum to l_{max}=3"

