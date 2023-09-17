#!/usr/bin/env gnuplot
dt = 0.0883795088135536 # Wce
dx = 0.5 # de
B0 = 0.0008018755454710313

set terminal pdfcairo enhanced color dashed font "Helvetica,10"
set output "figure1.pdf"

set xlabel "x (d_e)"
set ylabel "t {/Symbol W}_{ce}"
set cblabel "B_z / B_0"
set xrange [0:500]
set yrange [0:]
set cbrange [-0.1:0.1]
set palette defined (-1 "red", 0 "gray", 1 "blue")
set grid front
unset grid
set format cb "%.2f"
plot "run0_1d3v/B2.dat" using ($1*dx):($2*dt):($3/B0) matrix with image notitle
