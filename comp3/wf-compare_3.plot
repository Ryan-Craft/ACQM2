
set grid
set title "Numerov Method vs Analytical wf's in the Quantum Harmonic Oscillator Potential, n=3"

v(x) = 0.5*x**2

psi_3(x) = pi**(-0.25)*exp(-x**2 / 2)*(1/sqrt(3))*(2*x**3-3*x) 

set xrange[-5:5]
set yrange[2.5:4.5]

set xlabel "r (atomic units)"
set ylabel "Energy (hartrees)"

plot v(x) t "V(r)", psi_3(x)+3.5 t "Analytical psi_3(r)", "dx=0.5/psi_3.txt" u 1:($4+3.4833258359458590) with line t "psi_3(r) dr=0.5", "dx=0.1/psi_3.txt" u 1:($4+3.4999753371765463) with line t "psi_3(r) dr=0.1", "dx=0.01/psi_3.txt" u 1:($4+3.4999999986659764) with line t "psi_3(r) dr=0.01"

