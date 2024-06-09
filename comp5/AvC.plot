
set grid

set title font ",12"
set xlabel font ",12"
set ylabel font ",12"
set tics font ",12"
set key font ",12"

set title "10 eV Analytical vs Calculated Direct and\n Exchange on-shell matrix elements 1s-3s"

set ylabel "V(k',k)"
set xlabel "k"

plot "DirectAvsC.txt" u 1:2 t "Direct Analytical" w line, "DirectAvsC.txt" u 1:3 t "Calculated Direct" w line, "DirectAvsC.txt" u 1:4 t "Exchange Analytical"w line, "DirectAvsC.txt" u 1:5 t "Calculated Exchange" w line
