#!/usr/bin/gnuplot

set terminal pngcairo

set xlabel "{/Symbol w}_{pe} t"
set ylabel "{/Symbol D}E / E_0"

E1 = `sed -n -e 4p thermal_equib3_1/energy.dat | awk '{print $2}'`
E3 = `sed -n -e 4p thermal_equib3_3/energy.dat | awk '{print $2}'`
E5 = `sed -n -e 4p thermal_equib3_5/energy.dat | awk '{print $2}'`
E10 = `sed -n -e 4p thermal_equib3_10/energy.dat | awk '{print $2}'`
E30 = `sed -n -e 4p thermal_equib3_30/energy.dat | awk '{print $2}'`
E100 = `sed -n -e 4p thermal_equib3_100/energy.dat | awk '{print $2}'`
E300 = `sed -n -e 4p thermal_equib3_300/energy.dat | awk '{print $2}'`
E1000 = `sed -n -e 4p thermal_equib3_1000/energy.dat | awk '{print $2}'`

set title "nppc = 50"
set output "E_50.png"
plot "thermal_equib3_1/energy.dat"    u 1:($2/E1-1.)    w l t "Xi 1"   \
,    "thermal_equib3_3/energy.dat"    u 1:($2/E3-1.)    w l t "Xi 3"   \
,    "thermal_equib3_5/energy.dat"    u 1:($2/E5-1.)    w l t "Xi 5"   \
,    "thermal_equib3_10/energy.dat"   u 1:($2/E10-1.)   w l t "Xi 10"  \
,    "thermal_equib3_30/energy.dat"   u 1:($2/E30-1.)   w l t "Xi 30"  \
,    "thermal_equib3_100/energy.dat"  u 1:($2/E100-1.)  w l t "Xi 100" \
,    "thermal_equib3_300/energy.dat"  u 1:($2/E300-1.)  w l t "Xi 300" \
,    "thermal_equib3_1000/energy.dat" u 1:($2/E1000-1.) w l t "Xi 1000"


set title "nppc = 500"
set output "E_500.png"
E1 = `sed -n -e 4p thermal_equib_1/energy.dat | awk '{print $2}'`
E3 = `sed -n -e 4p thermal_equib_3/energy.dat | awk '{print $2}'`
E10 = `sed -n -e 4p thermal_equib_10/energy.dat | awk '{print $2}'`
E30 = `sed -n -e 4p thermal_equib_30/energy.dat | awk '{print $2}'`
E100 = `sed -n -e 4p thermal_equib_100/energy.dat | awk '{print $2}'`
E300 = `sed -n -e 4p thermal_equib_300/energy.dat | awk '{print $2}'`
E1000 = `sed -n -e 4p thermal_equib_1000/energy.dat | awk '{print $2}'`

plot "thermal_equib_1/energy.dat"    u 1:($2/E1-1.)    w l t "Xi 1"   \
,    "thermal_equib_3/energy.dat"    u 1:($2/E3-1.)    w l t "Xi 3"   \
,    "thermal_equib_10/energy.dat"   u 1:($2/E10-1.)   w l t "Xi 10"  \
,    "thermal_equib_30/energy.dat"   u 1:($2/E30-1.)   w l t "Xi 30"  \
,    "thermal_equib_100/energy.dat"  u 1:($2/E100-1.)  w l t "Xi 100" \
,    "thermal_equib_300/energy.dat"  u 1:($2/E300-1.)  w l t "Xi 300" \
,    "thermal_equib_1000/energy.dat" u 1:($2/E1000-1.) w l t "Xi 1000"


set title "nppc = 5000"
set output "E_5000.png"
E1 = `sed -n -e 4p thermal_equib2_1/energy.dat | awk '{print $2}'`
E3 = `sed -n -e 4p thermal_equib2_3/energy.dat | awk '{print $2}'`
E5 = `sed -n -e 4p thermal_equib2_5/energy.dat | awk '{print $2}'`
E10 = `sed -n -e 4p thermal_equib2_10/energy.dat | awk '{print $2}'`
E30 = `sed -n -e 4p thermal_equib2_30/energy.dat | awk '{print $2}'`
E100 = `sed -n -e 4p thermal_equib2_100/energy.dat | awk '{print $2}'`
E300 = `sed -n -e 4p thermal_equib2_300/energy.dat | awk '{print $2}'`
E1000 = `sed -n -e 4p thermal_equib2_1000/energy.dat | awk '{print $2}'`

plot "thermal_equib2_1/energy.dat"    u 1:($2/E1-1.)    w l t "Xi 1"   \
,    "thermal_equib2_3/energy.dat"    u 1:($2/E3-1.)    w l t "Xi 3"   \
,    "thermal_equib2_5/energy.dat"    u 1:($2/E5-1.)    w l t "Xi 5"   \
,    "thermal_equib2_10/energy.dat"   u 1:($2/E10-1.)   w l t "Xi 10"  \
,    "thermal_equib2_30/energy.dat"   u 1:($2/E30-1.)   w l t "Xi 30"  \
,    "thermal_equib2_100/energy.dat"  u 1:($2/E100-1.)  w l t "Xi 100" \
,    "thermal_equib2_300/energy.dat"  u 1:($2/E300-1.)  w l t "Xi 300" \
,    "thermal_equib2_1000/energy.dat" u 1:($2/E1000-1.) w l t "Xi 1000"

