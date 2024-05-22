

set title "Calculated 1ssg/2psu State vs Accurate Values from PEC Files\n N=10, alpha=1, lmax = 4"
set title font ",12"
set xlabel "Internuclear Separation (a_0)"
set xlabel font ",12"
set ylabel "Energy (Hartrees)"
set ylabel font ",12"

set tics font ",12"
set grid
set yrange[-1:4]
set xrange[0:8]
set pointsize 1

plot 'PEC.1' u 1:3 t "Calculated 1ssg" pointtype 7, 'PEC.2' u 1:3 t "Calculated 2psu" pointtype 7,'PEC_good/PEC.1ssg' t "PEC.1ssg" w line linecolor 1, 'PEC_good/PEC.2psu' t "PEC.2psu" w line linecolor 2

