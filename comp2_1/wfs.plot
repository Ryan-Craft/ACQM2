

set yrange[-0.61:-0.54]
set xrange[0:5]
set key font ",12"
set key top right

scale=0.005

set title "First Six H_2^+ Isotope Wavefunctions Inside the 1ssg Potential"
set title font ",12"
set xlabel "Internuclear Separation (a_0)"
set xlabel font ",12"
set ylabel "Energy (Hartrees)"
set ylabel font ",12"
set tics font ",11"


plot 'PEC_good/PEC.1ssg' u 1:2 w line t "1ssg electron potential", 'wfout.txt' u 1:(scale*$2-0.600) w line t "ψ_1", 'wfout.txt' u 1:(scale*$3-0.594) w line t "ψ_2", 'wfout.txt' u 1:(scale*$4-0.588) w line t "ψ_3", 'wfout.txt' u 1:(scale*$5-0.582) w line t "ψ_4" 

