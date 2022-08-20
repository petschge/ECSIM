#!/usr/bin/gnuplot

set terminal pdfcairo

set xrange [0:125]
set yrange [1e-1:1e5]
set xlabel "{/Symbol w}_{pe} t"
set ylabel "Electrical energy"

set format y "10^{%T}"
set logscale y

set output "1d1v_noB/figure1.pdf"
plot "1d1v_noB/energy.dat" u 1:($10+$11+$12) w l t "ECIM" \
,    10*exp(2*0.35*(x-10))

set output "1d1v/figure1.pdf"
plot "1d1v/energy.dat" u 1:($10+$11+$12) w l t "ECIM" \
,    10*exp(2*0.35*(x-10))

set output "1d3v/figure1.pdf"
plot "1d3v/energy.dat" u 1:($10+$11+$12) w l t "ECIM" \
,    10*exp(2*0.35*(x-10))




set output "figure1_comp.pdf"
plot "1d1v_noB/energy.dat" u 1:($10+$11+$12) w l t "ECIM, 1d1v, no B" \
,    "1d1v/energy.dat" u 1:($10+$11+$12) w l t "ECIM, 1d1v" \
,    "1d3v/energy.dat" u 1:($10+$11+$12) w l t "ECIM, 1d3v" \
,    10*exp(2*0.35*(x-10)) w l t "theory"
