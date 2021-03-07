
TPIBA=1.1845361

set pm3d map interpolate 2,2
set yrange [6:12]
splot './SC/BaTiO_Unfolded_SPFN.dat' binary record=( 121,-1) format='%double' u ($1/TPIBA):2:3 notitle, './BC/BaTiO3.bands.dat' u 1:2:2 w p ps 0.1 notitle
