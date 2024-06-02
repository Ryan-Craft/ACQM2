
set grid

set xrange [0:2]
set yrange [-1.5:1.5]

plot "threefuns.txt" u 1:2 w line t "potential V", "threefuns.txt" u 1:3 w line t "contwaves(2,:)", "threefuns.txt" u 1:4 w line t "contwaves(3,:)" 



