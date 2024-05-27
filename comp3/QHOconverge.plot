
set grid
set title "Numerov Method vs Analytical wf's in the Quantum Harmonic Oscillator Potential"

v(x) = 0.5*x**2

psi_0(x) = pi**(-0.25)*exp(-x**2 / 2) 
psi_1(x) = pi**(-0.25)* exp(-x**2 / 2) * sqrt(2)*x
psi_2(x) = pi**(-0.25)*exp(-x**2 / 2) * (1/sqrt(2))*(2*x**2 - 1)
psi_3(x) = pi**(-0.25)*exp(-x**2 / 2) * (1/sqrt(3)) * (2*x**3 - 3*x)


set xrange[-10:10]
set yrange[-2:5]

plot v(x), "psi_0.txt" u 1:4 with line, "psi_1.txt" u 1:4 with line, "psi_2.txt" u 1:4 with line, "psi_3.txt" u 1:4 with line

