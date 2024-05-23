
set grid
set title "Numerov Method vs Analytical wf's in the Quantum Harmonic Oscillator Potential"

v(x) = 0.5*x**2

psi_0(x) = pi**(-0.25)*exp(-x**2 / 2) 
psi_1(x) = pi**(-0.25)* exp(-x**2 / 2) * sqrt(2)*x
psi_2(x) = pi**(-0.25)*exp(-x**2 / 2) * (1/sqrt(2))*(2*x**2 - 1)
psi_3(x) = pi**(-0.25)*exp(-x**2 / 2) * (1/sqrt(3)) * (2*x**3 - 3*x)


set xrange[-7:7]
set yrange[0:5]

plot v(x), psi_0(x)+0.5, "dx=0.5/psi_0.txt" u 1:($4+0.49974711028260788) with line, "dx=0.1/psi_0.txt" u 1:($4+0.49999968333074568) with line, "dx=0.01/psi_0.txt" u 1:($4+0.50000000419848212) with line, "dx=0.5/psi_1.txt" u 1:($4+1.4982089671345686) with line, "dx=0.5/psi_2.txt" u 1:($4+2.4935052012437571) with line, "dx=0.5/psi_3.txt" u 1:($4+3.4833258359458590) with line, "dx=0.1/psi_1.txt" u 1:($4+1.4999972647164823) with line, "dx=0.1/psi_2.txt" u 1:($4+2.4999901984807384) with line, "dx=0.1/psi_3.txt" u 1:($4+3.4999753371765463) with line, "dx=0.01/psi_1.txt" u 1:($4+1.4999999999057003) with line, "dx=0.01/psi_2.txt" u 1:($4+2.4999999967839908) with line, "dx=0.01/psi_3.txt" u 1:($4+3.4999999986659764) with line

