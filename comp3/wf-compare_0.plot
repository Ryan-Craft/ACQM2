
set grid
set title "Numerov Method vs Analytical wf's in the Quantum Harmonic Oscillator Potential, n=0"

v(x) = 0.5*x**2

psi_0(x) = pi**(-0.25)*exp(-x**2 / 2) 

set xlabel "r (atomic units)"
set ylabel "Energy (hartrees)"


set xrange[-4:4]
set yrange[0.25:1.3]

plot v(x) t "V(r)", psi_0(x)+0.5 t "Analytical psi_0(r)", "dx=0.5/psi_0.txt" u 1:($4+0.49974711028260788) with line t "Computed psi_0(r), dr=0.5", "dx=0.1/psi_0.txt" u 1:($4+0.49999968333074568) with line t "psi_0(r), dr=0.1", "dx=0.01/psi_0.txt" u 1:($4+0.50000000419848212) with line t "psi_0(r), dr=0.01"

