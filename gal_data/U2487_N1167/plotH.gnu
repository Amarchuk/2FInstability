#! /usr/bin/gnuplot -persist
unset title
unset key
set terminal postscript eps enhanced solid "Helvetica" 14
set output "finalHz.eps"
set xlabel "/pi R, arcsec"
set ylabel "Величина СКО, отн.ед" font "Helvetica,18"
set yrange [0:1]
set style line 1 lt 1 pt 7
plot "~/RMSresult" using 2 title "СКО" with linespoints linestyle 1 notitle
