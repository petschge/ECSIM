#!/usr/bin/gnuplot

set terminal pngcairo color enhanced size 920,480
set xlabel "T (1 / ω)"
set ylabel "E (a.u.)"
set logscale y
set format y "10^{%T}"
set key top rmargin

set output "energy.png"
plot '1d3v/energy.dat' using 1:2  with lines title "Total energy"\
,    '1d3v/energy.dat' using 1:5  with lines title "B fields"    \
,    '1d3v/energy.dat' using 1:6  with lines title "E field"     \
,    '1d3v/energy.dat' using 1:3  with lines title "particles"   \
,    '1d3v/energy.dat' using 1:7  with lines title "B_x"         \
,    '1d3v/energy.dat' using 1:8  with lines title "B_y"         \
,    '1d3v/energy.dat' using 1:9  with lines title "B_z"         \
,    '1d3v/energy.dat' using 1:10 with lines title "E_x"         \
,    '1d3v/energy.dat' using 1:11 with lines title "E_y"         \
,    '1d3v/energy.dat' using 1:12 with lines title "E_z"

