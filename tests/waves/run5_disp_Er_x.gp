#!/usr/bin/env gnuplot
set output 'run5_disp_Er_x.pdf'
set terminal pdfcairo enhanced color dashed font "Helvetica,10"
set logscale cb
set cblabel "E (erg cm^{-3})" offset 1,0
set format cb "10^{%L}"
set tics front
set title "E_r in x direction (transverse)"
set ylabel '{/Symbol w} ({/Symbol W}_{i})'
set xlabel 'k (1/r_i)'

n_output_timesteps = 3535
n_timesteps = 3534
n_lineout_raumschritte = 1
nx = 3415

dt = 1.1316198311594571e-09 #s
dx = 4.798175416133162 #cm
c = 29979000000.0 #cm/s
omegape = 1000000000.0 #rad/s
omegapi = 233380014.0046683 #rad/s
Omegae = 500003584.6698209 #rad/s
Omegai = 27233310.711863883 #rad/s

k_nyquist = 0.6545560743287713 #1/cm
kalias_em = 0.08641633261881844 #1/cm
normk = 0.05707726666354434 #1/cm
normfreq = 27233310.711863883 #rad/s

kprime(k) = asin(k * dx/2.)*2./dx #Finite grid effects
oprime(o) = 2./dt * asin(dt/2. * o) #Finite grid effects
#kprime(k) = k
#oprime(o) = o

set yrange [0:1.5]
set xrange [-0.5:0.5]
set cbrange [1e-1:1e4]
set samples nx

plot "run5_1d3v/disp_Er_x.dat" matrix using (($1-nx/2)*2*pi/dx/nx/normk/n_lineout_raumschritte):($2*2*pi/dt/normfreq/n_timesteps):($3*n_lineout_raumschritte*n_lineout_raumschritte) with image notitle \
, 'run5_lmode' using (kprime($1)/normk):(oprime($2)/normfreq) w l lt 1 lc rgb "red" title "L mode" \
, 'run5_rmode' using (kprime($1)/normk):(oprime($2)/normfreq) w l lt 1 lc rgb "green" title "R mode" \
, sqrt(omegape**2 + omegapi**2)/normfreq lt 2 lc rgb "black" title "plasma freq {/Symbol w}_{p}" \
, Omegae/normfreq lt 2 lc rgb "green" title "electron gyro freq {/Symbol W}_{e}" \
, Omegai/normfreq lt 2 lc rgb "green" title "ion gyro freq {/Symbol W}_{i}" \
