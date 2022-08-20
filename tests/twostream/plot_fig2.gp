#!/usr/bin/gnuplot

set terminal pdfcairo

set xrange [0:125]
set yrange [1e-16:1e0]
set xlabel "{/Symbol w}_{pe} t"
set ylabel "error in total energy: |E(t) - E_0| / E_0"

set format y "10^{%T}"
set logscale y

set output "1d1v_noB/figure2.pdf"
E0 = `sed -n -e 3p 1d1v_noB/energy.dat | awk '{print $2}'`
plot "1d1v_noB/energy.dat" u 1:(sqrt(($2-E0)**2)/E0) w l t "ECIM"

set output "1d1v/figure2.pdf"
E1 = `sed -n -e 4p 1d1v/energy.dat | awk '{print $2}'`
plot "1d1v/energy.dat" u 1:(sqrt(($2-E1)**2)/E1) w l t "ECIM"

set output "1d3v/figure2.pdf"
E2 = `sed -n -e 4p 1d3v/energy.dat | awk '{print $2}'`
plot "1d3v/energy.dat" u 1:(sqrt(($2-E2)**2)/E2) w l t "ECIM"





set output "figure2.pdf"
plot "1d1v_noB/energy.dat" u 1:(sqrt(($2-E0)**2)/E0) w l t "ECIM, 1d1v, no B" \
,    "1d1v/energy.dat" u 1:(sqrt(($2-E1)**2)/E1) w l t "ECIM, 1d1v" \
,    "1d3v/energy.dat" u 1:(sqrt(($2-E2)**2)/E2) w l t "ECIM, 1d3v"
