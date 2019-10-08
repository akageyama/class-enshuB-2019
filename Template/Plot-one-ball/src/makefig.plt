set size square
set title "Particle"
set xrange [0:1.0]
set yrange [0:1.0]
set nokey

plot 'particle.pos.data' with p ps 3 pt 15 lc 'blue'

replot
pause -1

