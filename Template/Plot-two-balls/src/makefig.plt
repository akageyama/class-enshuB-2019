set size square
set title "Particles"
set xrange [0:1.0]
set yrange [0:1.0]
set nokey

plot "particles.pos.data" with p ps 3 pt 15 lc 'blue'

pause -1

