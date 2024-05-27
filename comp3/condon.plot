set grid

set title font ",12"
set tics font ",11"
set key font ",12"
set xlabel "Energy (eV)" font ",12"
set ylabel "KER Distribution (abritrary units)" font ",12"

set title "H_2^+ 1ssg -> 2psu Transition KER Distribution Using Frank-Condon Approximation"

set xtics "2"

set xrange [0:25]


plot "condonout.txt" u 1:2 with line t "v=0", "condonout.txt" u 1:3 with line t "v=3","condonout.txt" u 1:4 with line t "v=6","condonout.txt" u 1:5 with line lc 7 t "v=9"


