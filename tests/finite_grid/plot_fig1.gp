#!/usr/bin/gnuplot

set terminal pdfcairo

set logscale x
set logscale y
set format x "10^{%T}"
set format y "10^{%T}"
set key top right
set yrange [1e-30:1]
set xrange [1e1:1e16]

set xlabel "{/Symbol D}x / {/Symbol l}_D"

set output "figure1_1d1v_noB.pdf"
set title "1d1v, no B, {/Symbol w}_{pe} {/Symbol D}t = 0.0125 {/Symbol X}"
plot "Xi_1d1v_noB.dat" u ($1*10):2             w l  t "v_{th}/c"\
,    ""       u ($1*10):(sqrt($3*$3)) w l  t "<{/Symbol D}E_t>/E_t(0)" \
,    ""       u ($1*10):4             w l  t "{/Symbol d}E_E / E_t(0)" \
,    ""       u ($1*10):5             w l  t "{/Symbol d}E_K / E_t(0)" \
,    ""       u ($1*10):6             w l  t "{/Symbol n}/{/Symbol w}_{pe}"

set output "figure1_1d1v.pdf"
set title "1d1v, {/Symbol w}_{pe} {/Symbol D}t = 0.0125 {/Symbol X}"
plot "Xi_1d1v.dat" u ($1*10):2             w l  t "v_{th}/c"\
,    ""       u ($1*10):(sqrt($3*$3)) w l  t "<{/Symbol D}E_t>/E_t(0)" \
,    ""       u ($1*10):4             w l  t "{/Symbol d}E_E / E_t(0)" \
,    ""       u ($1*10):5             w l  t "{/Symbol d}E_K / E_t(0)" \
,    ""       u ($1*10):6             w l  t "{/Symbol n}/{/Symbol w}_{pe}"

set output "figure1_1d3v.pdf"
set title "1d3v, {/Symbol w}_{pe} {/Symbol D}t = 0.0125 {/Symbol X}"
plot "Xi_1d3v.dat" u ($1*10):2             w l  t "v_{th}/c"\
,    ""       u ($1*10):(sqrt($3*$3)) w l  t "<{/Symbol D}E_t>/E_t(0)" \
,    ""       u ($1*10):4             w l  t "{/Symbol d}E_E / E_t(0)" \
,    ""       u ($1*10):5             w l  t "{/Symbol d}E_K / E_t(0)" \
,    ""       u ($1*10):6             w l  t "{/Symbol n}/{/Symbol w}_{pe}"
