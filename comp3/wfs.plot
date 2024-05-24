

set yrange[-0.6:1]
set xrange[0:10]
set key font ",12"
set key top right

scale=0.05

set title "E_k = 5, 10, 15, 20 eV Dissociated H_2^+ Wavefunction"
set title font ",12"
set xlabel "Internuclear Separation (a_0)"
set xlabel font ",12"
set ylabel "Energy (Hartrees)"
set ylabel font ",12"
set tics font ",11"


plot 'PEC_good/PEC.2psu' u 1:2 w line, "unboundout.txt" u 1:(scale*$3-0.316) w line t "Ek=5eV", "unboundout.txt" u 1:(scale*$4-0.1325) w line t "Ek=10eV", "unboundout.txt" u 1:(scale*$5+0.0512) w line t "Ek=15eV", "unboundout.txt" u 1:(scale*$5+0.235) w line t "Ek=20eV"

