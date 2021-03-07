set yrange [-5:1]
set arrow from  0.422104,   -5 to  0.422104,  1 nohead
set arrow from  0.844208,   -5 to  0.844208,  1 nohead
set arrow from  1.441153,   -5 to  1.441153,  1 nohead
set arrow from  2.038098,   -5 to  2.038098,  1 nohead
set arrow from  2.460202,   -5 to  2.460202,  1 nohead
set arrow from  3.057147,   -5 to  3.057147,  1 nohead
set arrow from  3.479251,   -5 to  3.479251,  1 nohead
set xtics (" G "  0.000000, " X "  0.422104, " R "  0.844208, " G "  1.441153, " M "  2.038098, " A "  2.460202, " Z "  3.057147, " G "  3.479251)

EFERMI=11.071896
set ylabel "Energy [eV]" textcolor rgb text_color font my_font
plot 'BaTiO3.bands.dat' u 1:($2-EFERMI):(1+$2*0) w p notitle ps 0.5 pt 6 palette
