
set grid
set title "Numerov Method vs Analytical wf's in the Quantum Harmonic Oscillator Potential, n=1"

v(x) = 0.5*x**2

psi_1(x) = pi**(-0.25)*exp(-x**2 / 2) *sqrt(2)*x

set xrange[-4:4]
set yrange[0.5:3.2]

set xlabel "r (atomic units)"
set ylabel "Energy (hartrees)"

plot v(x) t "V(r)", psi_1(x)+1.5 t "Analytical psi_1(r)", "dx=0.5/psi_1.txt" u 1:($4+1.4982089671345686) with line t "Computed psi_1(r), dr=0.5", "dx=0.1/psi_1.txt" u 1:($4+1.4999972647164823) with line t "psi_1(r), dr=0.1", "dx=0.01/psi_1.txt" u 1:($4+1.4999999999057003) with line t "psi_1(r), dr=0.01"

