#!/usr/bin/env gnuplot
dt = 0.0883795088135536 # Wce
dx = 0.5 # de
Nx = 1000
Nt = 10000

kmax = pi/dx
omegamax = pi/dt


set terminal pdfcairo enhanced color dashed font "Helvetica,10"
set output "figure2.pdf"

set xlabel "k_x d_e"
set ylabel "{/Symbol w} / {/Symbol W}_{ce}"
set xrange [-3:3]
set yrange [0:1.5]
set cbrange [5e-8:5.]
set grid front
unset grid
set format y "%.1f"
set format cb "10^{%T}"

set logscale cb
plot "run0_1d3v/disp_Er_x.dat" using (2.*(0.5-$1/Nx)*kmax):(2.*$2/Nt*omegamax):3 matrix with image notitle
