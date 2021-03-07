
TPIBA=1.159259

set pm3d map interpolate 2,2
set yrange [10:20]
splot './SC/Fe-unfold.dat1' binary record=(2001,-1) format='%double' u ($1*TPIBA):2:3 notitle, './BC/Fe.bnd_UP.dat' u 1:2:2 w p ps 0.1 notitle
