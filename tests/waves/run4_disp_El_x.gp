#!/usr/bin/env gnuplot
set output 'run4_disp_El_x.pdf'
set terminal pdfcairo enhanced color dashed font "Helvetica,10"
set logscale cb
set cblabel "E (erg cm^{-3})" offset 1,0
set format cb "10^{%L}"
set tics front
set title "E_l in x direction (transverse)"
set ylabel '{/Symbol w} ({/Symbol W}_{e})'
set xlabel 'k ({/Symbol p}/d_e)'

n_output_timesteps = 2561
n_timesteps = 2560
n_lineout_raumschritte = 1
nx = 512

dt = 3.7735004930567603e-10 #s
dx = 16.0 #cm
c = 29979000000.0 #cm/s
omegape = 1000000000.0 #rad/s
omegapi = 23338001.40046683 #rad/s
Omegae = 500035243.60188156 #rad/s
Omegai = 272350.35054568714 #rad/s

k_nyquist = 0.19596604565239065 #1/cm
kalias_em = 0.275806640156716 #1/cm
normk = 0.10479311029686758 #1/cm
normfreq = 500035243.60188156 #rad/s

kprime(k) = asin(k * dx/2.)*2./dx #Finite grid effects
oprime(o) = 2./dt * asin(dt/2. * o) #Finite grid effects
#kprime(k) = k
#oprime(o) = o

set yrange [0:1.5]
set xrange [-1.0:1.0]
set cbrange [1e-1:1e4]
set samples nx

plot "run4_1d3v/disp_El_x.dat" matrix using (($1-nx/2)*2*pi/dx/nx/normk/n_lineout_raumschritte):($2*2*pi/dt/normfreq/n_timesteps):($3*n_lineout_raumschritte*n_lineout_raumschritte) with image notitle \
, 'lmode' using (kprime($1)/normk):(oprime($2)/normfreq) w l lt 1 lc rgb "red" title "L mode" \
, 'rmode' using (kprime($1)/normk):(oprime($2)/normfreq) w l lt 1 lc rgb "green" title "R mode" \
, 'xmode' using (kprime($1)/normk):(oprime($2)/normfreq) w l lt 1 lc rgb "blue" title "X mode" \
, sqrt(omegape**2 + omegapi**2)/normfreq lt 2 lc rgb "black" title "plasma freq {/Symbol w}_{p}" \
, Omegae/normfreq lt 2 lc rgb "green" title "electron gyro freq {/Symbol W}_{e}" \
, Omegai/normfreq lt 2 lc rgb "green" title "ion gyro freq {/Symbol W}_{i}" \
