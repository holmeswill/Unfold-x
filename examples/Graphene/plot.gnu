


#set term png enhanced large
#set output plot.gnu
EF=-1.7106
tpiba=0.5 #output from bands.x is already in tpiba (so 1/2 is BCalat/SCalat)

set pm3d map interpolate 2,2

set xtics ("K" 0.0, "{/Symbol G}" 0.666666, "M" 1.244, "K" 1.577349)
set ylabel "E - E_F [eV]"

unset key


splot [:][-20:8]'./SC/graphene.dat' binary record=(3001,-1) format='%double' u ($1*tpiba):($2-EF):3,'./BC/graphene.bands.dat' u 1:($2-EF):(4) w p ps 0.5
