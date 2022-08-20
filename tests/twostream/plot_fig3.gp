#!/usr/bin/gnuplot

set terminal pdfcairo size 5.,2.

set xrange [0:125]
set xlabel "{/Symbol w}_{pe} t"

#set format y "10^{%T}"
#set logscale y
set format y "%.3f"
set ytics 1e-3
set mytics 2

# derivative functions.  Return 1/0 for first point, otherwise delta y or (delta y)/(delta x)
d(y) = ($0 == 0) ? (y1 = y, 1/0) : (y2 = y1, y1 = y, y1-y2)
d2(x,y) = ($0 == 0) ? (x1 = x, y1 = y, 1/0) : (x2 = x1, x1 = x, y2 = y1, y1 = y, (y1-y2)/(x1-x2))

set output "1d1v_noB/figure3.pdf"
set ylabel " P_x(t) / P_{x,0}"
set yrange [-3e-3:3e-3]
plot "1d1v_noB/momentum.dat" u 1:($3/$2) w l t "ECIM, 1d1v, no B"

print "ECIM, 1d1v_noB"
stats [0:125] "1d1v_noB/momentum.dat" u 1:(d2($1,$3)/$2) name "all_noB" nooutput
print "nu_rms during entire time ",all_noB_stddev_y," w_pe"
stats [20:125] "1d1v_noB/momentum.dat" u 1:(d2($1,$3)/$2) name "instability_noB" nooutput
print "nu_rms during instability ",instability_noB_stddev_y," w_pe"

set output "1d1v_noB/figure3b.pdf"
set ylabel "{/Symbol n} / {/Symbol w}_{pe}"
set yrange [-4e-3:4e-3]
plot "1d1v_noB/momentum.dat" u 1:(d2($1,$3)/$2) w l t sprintf("ECIM, 1d1v, no B: <{/Symbol n}> = %.2e {/Symbol w}_{pe}", instability_noB_stddev_y)






set output "1d1v/figure3.pdf"
set ylabel " P_x(t) / P_{x,0}"
set yrange [-3e-3:3e-3]
plot "1d1v/momentum.dat" u 1:($3/$2) w l t "ECIM, 1d1v"

print "ECIM, 1d1v"
stats [0:125] "1d1v/momentum.dat" u 1:(d2($1,$3)/$2) name "all_1d1v" nooutput
print "nu_rms during entire time ",all_1d1v_stddev_y," w_pe"
stats [20:125] "1d1v/momentum.dat" u 1:(d2($1,$3)/$2) name "instability_1d1v" nooutput
print "nu_rms during instability ",instability_1d1v_stddev_y," w_pe"

set output "1d1v/figure3b.pdf"
set ylabel "{/Symbol n} / {/Symbol w}_{pe}"
set yrange [-4e-3:4e-3]
plot "1d1v/momentum.dat" u 1:(d2($1,$3)/$2) w l t sprintf("ECIM, 1d1v: <{/Symbol n}> = %.2e {/Symbol w}_{pe}", instability_1d1v_stddev_y)





set output "1d3v/figure3.pdf"
set ylabel " P_x(t) / P_{x,0}"
set yrange [-3e-3:3e-3]
plot "1d3v/momentum.dat" u 1:($3/$2) w l t "ECIM, 1d3v"

print "ECIM, 1d3v"
stats [0:125] "1d3v/momentum.dat" u 1:(d2($1,$3)/$2) name "all_1d3v" nooutput
print "nu_rms during entire time ",all_1d3v_stddev_y," w_pe"
stats [20:125] "1d3v/momentum.dat" u 1:(d2($1,$3)/$2) name "instability_1d3v" nooutput
print "nu_rms during instability ",instability_1d3v_stddev_y," w_pe"

set output "1d3v/figure3b.pdf"
set ylabel "{/Symbol n} / {/Symbol w}_{pe}"
set yrange [-4e-3:4e-3]
plot "1d3v/momentum.dat" u 1:(d2($1,$3)/$2) w l t sprintf("ECIM, 1d3v: <{/Symbol n}> = %.2e {/Symbol w}_{pe}", instability_1d3v_stddev_y)






set output "figure3.pdf"
set ylabel " P_x(t) / P_{x,0}"
set yrange [-3e-3:3e-3]
plot "1d1v_noB/momentum.dat" u 1:($3/$2) w l t "ECIM, 1d1v, no B" \
,    "1d1v/momentum.dat"     u 1:($3/$2) w l t "ECIM, 1d31" \
,    "1d3v/momentum.dat"     u 1:($3/$2) w l t "ECIM, 1d3v"

set output "figure3b.pdf"
set ylabel "{/Symbol n} / {/Symbol w}_{pe}"
set yrange [-4e-3:4e-3]
plot "1d1v_noB/momentum.dat" u 1:(d2($1,$3)/$2) w l t sprintf("ECIM, 1d1v, no B: <{/Symbol n}> = %.2e {/Symbol w}_{pe}", instability_noB_stddev_y) \
,    "1d1v/momentum.dat"     u 1:(d2($1,$3)/$2) w l t sprintf("ECIM, 1d1v: <{/Symbol n}> = %.2e {/Symbol w}_{pe}", instability_1d1v_stddev_y) \
,    "1d3v/momentum.dat"     u 1:(d2($1,$3)/$2) w l t sprintf("ECIM, 1d3v: <{/Symbol n}> = %.2e {/Symbol w}_{pe}", instability_1d3v_stddev_y)





