set size square
set title "Particles"
set xrange [0:1.0]
set yrange [0:1.0]
set nokey

max_frame = 100

do for [i=0: max_frame] {
  plot sprintf("particles.pos.data.%05.0f", i) with p ps 3 pt 15 lc 'blue'
  pause 0.1
}
pause -1

