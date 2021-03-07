set arrow from  0.730603,   -5 to  0.730603,  5 nohead
set arrow from  1.461206,   -5 to  1.461206,  5 nohead
set arrow from  2.494434,   -5 to  2.494434,  5 nohead
set arrow from  3.759876,   -5 to  3.759876,  5 nohead
set arrow from  4.793104,   -5 to  4.793104,  5 nohead
set arrow from  5.523707,   -5 to  5.523707,  5 nohead
set arrow from  6.254310,   -5 to  6.254310,  5 nohead
set xtics (" G "  0.000000, " X "  0.730603, " M "  1.461206, " G "  2.494434, " R "  3.759876, " X "  4.793104, " M "  5.523707, " R "  6.254310)

EFERMI=17.5042
plot 'Fe.bands.dat' u 1:($2-EFERMI) w l notitle
