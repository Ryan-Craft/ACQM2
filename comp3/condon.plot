set grid

set title font ",12"
set tics font ",11"

set xlabel "Energy (eV)" font ",12"
set ylabel "KER Distribution (abritrary units)" font ",12"

plot "condonout.txt" u 1:2 with line, "condonout.txt" u 1:3 with line ,"condonout.txt" u 1:4 with line ,"condonout.txt" u 1:2 with line ,"condonout.txt" u 1:5 with line


