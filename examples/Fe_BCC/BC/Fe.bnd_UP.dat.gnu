set arrow from  1.003948,   -50 to  1.003948,  25 nohead
set arrow from  1.823668,   -50 to  1.823668,  25 nohead
set arrow from  2.325642,   -50 to  2.325642,  25 nohead
set xtics (" G "  0.000000, " H "  1.003948, " N "  1.823668, " P "  2.325642)

EFERMI=0
plot 'Fe.bnb_UP.dat' u 1:($2-EFERMI) w l notitle
