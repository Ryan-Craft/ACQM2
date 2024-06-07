
set grid

splot 'Vdirectout.txt' u 1:2:3 with lines, 'Exchange.txt' u 1:2:3 with lines
set xyplane 0.1

set xrange [0:5]
set yrange [0:5]

