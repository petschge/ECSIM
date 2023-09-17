#!/usr/bin/env gnuplot

set terminal pdfcairo enhanced color dashed font "Helvetica,10"
set output "figure3.pdf"

wpe_Wce = 4

set xlabel "t {/Symbol W}_{ce}"
set ylabel "E"
set logscale y
set format y "10^{%T}"
set yrange [1e6:1e15]

plot "run0_1d3v/energy.dat" u ($1/wpe_Wce):2  w l t "total"    \
,    ""                     u ($1/wpe_Wce):7  w l t "B_x"      \
,    ""                     u ($1/wpe_Wce):8  w l t "B_y"      \
,    ""                     u ($1/wpe_Wce):9  w l t "B_z"      \
,    ""                     u ($1/wpe_Wce):10 w l t "E_x"      \
,    ""                     u ($1/wpe_Wce):11 w l t "E_y"      \
,    ""                     u ($1/wpe_Wce):12 w l t "E_z"      \
,    ""                     u ($1/wpe_Wce):13 w l t "electrons"\
,    ""                     u ($1/wpe_Wce):14 w l t "protons"

