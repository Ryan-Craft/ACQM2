
set grid
set title "Numerov Method vs Analytical wf's in the Quantum Harmonic Oscillator Potential, n=2"

v(x) = 0.5*x**2

psi_2(x) = pi**(-0.25)*exp(-x**2 / 2)*(1/sqrt(2))*(2*x**2-1) 

set xlabel "r (atomic units)"
set ylabel "Energy (hartrees)"

set xrange[-4:4]
set yrange[1.5:3.5]

plot v(x) t "V(r)", psi_2(x)+2.5 t "Analytical psi_2(r)", "dx=0.5/psi_2.txt" u 1:($4+2.4935052012437571) with line t "Computed psi_2(r), dr=0.5", "dx=0.1/psi_2.txt" u 1:($4+2.4999901984807384) with line t "psi_2(r), dr=0.1", "dx=0.01/psi_2.txt" u 1:($4+2.4999999967839908) with line t "psi_2(r), dr=0.01"

