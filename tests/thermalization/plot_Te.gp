#!/usr/bin/gnuplot

set terminal pngcairo

set xlabel "{/Symbol w}_{pe} t"
set ylabel "T_e / T_{e,0}"

E1 = `sed -n -e 4p thermal_equib3_1/energy.dat | awk '{print $13}'`
E3 = `sed -n -e 4p thermal_equib3_3/energy.dat | awk '{print $13}'`
E5 = `sed -n -e 4p thermal_equib3_5/energy.dat | awk '{print $13}'`
E10 = `sed -n -e 4p thermal_equib3_10/energy.dat | awk '{print $13}'`
E30 = `sed -n -e 4p thermal_equib3_30/energy.dat | awk '{print $13}'`
E100 = `sed -n -e 4p thermal_equib3_100/energy.dat | awk '{print $13}'`
E300 = `sed -n -e 4p thermal_equib3_300/energy.dat | awk '{print $13}'`
E1000 = `sed -n -e 4p thermal_equib3_1000/energy.dat | awk '{print $13}'`
Evpic50 = `sed -n -e 4p vpic/rundata_50/energy.dat | awk '{print $9}'`
dtvpic50 = `sed -n -e 3p vpic/rundata_50/energy.dat | cut -d = -f 2`
Evpic10000 = `sed -n -e 4p vpic/rundata_10000/energy.dat | awk '{print $9}'`
dtvpic10000 = `sed -n -e 3p vpic/rundata_10000/energy.dat | cut -d = -f 2`
Evpic410000 = `sed -n -e 4p vpic/rundata4_10000/energy.dat | awk '{print $9}'`
dtvpic410000 = `sed -n -e 3p vpic/rundata4_10000/energy.dat | cut -d = -f 2`

set title "nppc = 50"
set output "nppc_50.png"
set xrange [:3600.]
set yrange [0.5:5.]
plot "thermal_equib3_1/energy.dat"    u 1:($13/E1)    w l t "Xi 1"   \
,    "thermal_equib3_3/energy.dat"    u 1:($13/E3)    w l t "Xi 3"   \
,    "thermal_equib3_5/energy.dat"    u 1:($13/E5)    w l t "Xi 5"   \
,    "thermal_equib3_10/energy.dat"   u 1:($13/E10)   w l t "Xi 10"  \
,    "thermal_equib3_30/energy.dat"   u 1:($13/E30)   w l t "Xi 30"  \
,    "thermal_equib3_100/energy.dat"  u 1:($13/E100)  w l t "Xi 100" \
,    "thermal_equib3_300/energy.dat"  u 1:($13/E300)  w l t "Xi 300" \
,    "thermal_equib3_1000/energy.dat" u 1:($13/E1000) w l t "Xi 1000" \
,    "vpic/rundata_50/energy.dat" u ($1*dtvpic50):($9/Evpic50) every 1000 w l lw 2 dt 3 t "vpic 50" \
,    "vpic/rundata_10000/energy.dat" u ($1*dtvpic10000):($9/Evpic10000) every 1000 w l lw 2 dt 3 t "vpic 10000" \
,    "vpic/rundata4_10000/energy.dat" u ($1*dtvpic410000):($9/Evpic410000) every 1000 w l lw 2 dt 3 t "vpic 10000, Xi 4"

unset yrange

set title "nppc = 500"
set output "nppc_500.png"
E1 = `sed -n -e 4p thermal_equib_1/energy.dat | awk '{print $13}'`
E3 = `sed -n -e 4p thermal_equib_3/energy.dat | awk '{print $13}'`
E10 = `sed -n -e 4p thermal_equib_10/energy.dat | awk '{print $13}'`
E30 = `sed -n -e 4p thermal_equib_30/energy.dat | awk '{print $13}'`
E100 = `sed -n -e 4p thermal_equib_100/energy.dat | awk '{print $13}'`
E300 = `sed -n -e 4p thermal_equib_300/energy.dat | awk '{print $13}'`
E1000 = `sed -n -e 4p thermal_equib_1000/energy.dat | awk '{print $13}'`
Evpic500 = `sed -n -e 4p vpic/rundata_500/energy.dat | awk '{print $9}'`
dtvpic500 = `sed -n -e 3p vpic/rundata_500/energy.dat | cut -d = -f 2`
Evpic100000 = `sed -n -e 4p vpic/rundata_100000/energy.dat | awk '{print $9}'`
dtvpic100000 = `sed -n -e 3p vpic/rundata_100000/energy.dat | cut -d = -f 2`

plot "thermal_equib_1/energy.dat"    u 1:($13/E1)    w l t "Xi 1"   \
,    "thermal_equib_3/energy.dat"    u 1:($13/E3)    w l t "Xi 3"   \
,    "thermal_equib_10/energy.dat"   u 1:($13/E10)   w l t "Xi 10"  \
,    "thermal_equib_30/energy.dat"   u 1:($13/E30)   w l t "Xi 30"  \
,    "thermal_equib_100/energy.dat"  u 1:($13/E100)  w l t "Xi 100" \
,    "thermal_equib_300/energy.dat"  u 1:($13/E300)  w l t "Xi 300" \
,    "thermal_equib_1000/energy.dat" u 1:($13/E1000) w l t "Xi 1000" \
,    "vpic/rundata_500/energy.dat" u ($1*dtvpic500):($9/Evpic500) every 1000 w l lw 2 dt 3 t "vpic 500" \
,    "vpic/rundata_100000/energy.dat" u ($1*dtvpic100000):($9/Evpic100000) every 1000 w l lw 2 dt 3 t "vpic 100000"


set title "nppc = 5000"
set output "nppc_5000.png"
E1 = `sed -n -e 4p thermal_equib2_1/energy.dat | awk '{print $13}'`
E3 = `sed -n -e 4p thermal_equib2_3/energy.dat | awk '{print $13}'`
E5 = `sed -n -e 4p thermal_equib2_5/energy.dat | awk '{print $13}'`
E10 = `sed -n -e 4p thermal_equib2_10/energy.dat | awk '{print $13}'`
E30 = `sed -n -e 4p thermal_equib2_30/energy.dat | awk '{print $13}'`
E100 = `sed -n -e 4p thermal_equib2_100/energy.dat | awk '{print $13}'`
E300 = `sed -n -e 4p thermal_equib2_300/energy.dat | awk '{print $13}'`
E1000 = `sed -n -e 4p thermal_equib2_1000/energy.dat | awk '{print $13}'`
Evpic5000 = `sed -n -e 4p vpic/rundata_5000/energy.dat | awk '{print $9}'`
dtvpic5000 = `sed -n -e 3p vpic/rundata_5000/energy.dat | cut -d = -f 2`
Evpic1000000 = `sed -n -e 4p vpic/rundata_1000000/energy.dat | awk '{print $9}'`
dtvpic1000000 = `sed -n -e 3p vpic/rundata_1000000/energy.dat | cut -d = -f 2`
set xrange [0:1060]

plot "thermal_equib2_1/energy.dat"    u 1:($13/E1)    w l t "Xi 1"   \
,    "thermal_equib2_3/energy.dat"    u 1:($13/E3)    w l t "Xi 3"   \
,    "thermal_equib2_5/energy.dat"    u 1:($13/E5)    w l t "Xi 5"   \
,    "thermal_equib2_10/energy.dat"   u 1:($13/E10)   w l t "Xi 10"  \
,    "thermal_equib2_30/energy.dat"   u 1:($13/E30)   w l t "Xi 30"  \
,    "thermal_equib2_100/energy.dat"  u 1:($13/E100)  w l t "Xi 100" \
,    "thermal_equib2_300/energy.dat"  u 1:($13/E300)  w l t "Xi 300" \
,    "thermal_equib2_1000/energy.dat" u 1:($13/E1000) w l t "Xi 1000" \
,    "vpic/rundata_5000/energy.dat" u ($1*dtvpic5000):($9/Evpic5000) every 1000 w l lw 2 dt 3 t "vpic 5000" \
,    "vpic/rundata_1000000/energy.dat" u ($1*dtvpic1000000):($9/Evpic1000000) every 1000 w l lw 2 dt 3 t "vpic 1000000"

