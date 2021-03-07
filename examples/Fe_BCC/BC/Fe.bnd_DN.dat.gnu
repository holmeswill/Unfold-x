set arrow from  1.003948,   -5 to  1.003948,  5 nohead
set arrow from  1.823668,   -5 to  1.823668,  5 nohead
set arrow from  2.325642,   -5 to  2.325642,  5 nohead
set xtics (" G "  0.000000, " H "  1.003948, " N "  1.823668, " P "  2.325642)

EFERMI=0
plot 'Fe.bnb_DN.dat' u 1:($2-EFERMI) w l notitle
