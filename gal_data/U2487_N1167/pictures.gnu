set term post enhanced color solid "Arial" 20
set output 'vel.ps'

set xrange [-10:210]
set yrange [0:400]
set xlabel 'R, arcsec' 
set ylabel 'V_r, km/s'
set xtics 0,50,210 #scale 1.5,0.8
set mxtics 2
set ytics 0,100,400 #scale 1.5,0.8
set mytics 5

plot \
'v_gas.dat' using ($1):($3) title 'HI' w p pt 2 ps 1 lw 2,\
'v_gas.dat' using ($1):($3):($4) notitle w errorbars  pt 2 ps 1 lw 2,\
'v_stars_ma.dat' using ($1):((4960-$2)/sin(38*pi/180)) title "Stars" w p pt 5 ps 1 lw 4,\
'v_stars_ma.dat' using ($1):((4960-$2)/sin(38*pi/180)):($3) notitle w errorbars pt 5 ps 1 lw 4
