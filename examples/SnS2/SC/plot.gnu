set term pngcairo
set output "sns-npool.png"

set size 0.75 , 1.0
set origin 0.,0.

#set title "SnS2-npool"
set ylabel '{E} (eV)'
set noxlabel
set xrange [0:6]
set yrange [-10:15]

set pm3d map interpolate 2,2
splot './SnS2.dat' binary record=( 401,-1) format='%double' u 1:2:3

