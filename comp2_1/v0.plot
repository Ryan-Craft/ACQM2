

set yrange[0:3]
set xrange[0:7.5]
set key font ",12"
set key top right

scale=1
set grid

set title "First Four v=0 Vibrational Nuclear Wavefunctions Inside the 1ssg Potential"
set title font ",12"
set xlabel "Internuclear Separation (a_0)"
set xlabel font ",12"
set ylabel "Energy (Hartrees)"
set ylabel font ",12"
set tics font ",11"


plot 'H2.txt' u ($1-1):(scale*$2) w line t "H_2^+ wf", 'HD.txt' u ($1):(scale*$2) w line t "HD^+ wf", 'HT.txt' u ($1+1):(scale*$2) w line t "HT^+ wf", 'D2.txt' u ($1+2):(scale*$2) w line t "D_2^+ wf", 'DT.txt' u ($1+3):(scale*$2) w line t "DT^+" lc 8, 'T2.txt' u ($1+4):(scale*$2) w line t "T_2^+" 

