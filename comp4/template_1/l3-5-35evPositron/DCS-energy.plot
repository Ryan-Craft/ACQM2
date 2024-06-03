

set grid

set title font ",12"
set tics font ",11"

set xlabel font ",12"
set ylabel font ",12"

set xlabel "Scattering Angle (degrees)"
set ylabel "Differential Cross Section (a_0^2/Sr)"

set title "Differential Cross Section of Incident Energies 5 to 25 eV\n vs Scattering Angle"

plot "DCSout5.txt" u 1:2 t "Energy 5 eV", "DCSout10.txt" u 1:2 t "Energy 10 eV", "DCSout15.txt" u 1:2 t "Energy 15 eV", "DCSout20.txt" u 1:2 t "Energy 20 eV", "DCSout25.txt" u 1:2 t "Energy 25 eV", "DCSout30.txt" u 1:2 t "Energy 30 eV","DCSout35.txt" u 1:2 t "Energy 35 eV"
