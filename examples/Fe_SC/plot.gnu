TPIBA=1.0332292

set pm3d map interpolate 2,2
set yrange [10:20]
splot './SC/Fe.dat' binary record=(2001,-1) format='%double' u ($1*TPIBA):2:3 notitle, './BC/Fe.bands.dat' u 1:2:2 w p ps 0.1 notitle
