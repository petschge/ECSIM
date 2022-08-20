#!/usr/bin/env gnuplot
set output 'run0_disp_E2_x.pdf'
set terminal pdfcairo enhanced color dashed font "Helvetica,10"
set logscale cb
set cblabel "E (erg cm^{-3})" offset 1,0
set format cb "10^{%L}"
set tics front
set title "E_z in x direction (transverse)"
set ylabel '{/Symbol w} ({/Symbol w}_{pe})'
set xlabel 'k (1/{/Symbol l}_D)'

n_output_timesteps = 10401
n_timesteps = 10400
n_lineout_raumschritte = 1
nx = 2048

dt = 1.9230769230769234e-11 #s
dx = 1.0 #cm
c = 29979000000.0 #cm/s
omegape = 1000000000.0 #rad/s
omegapi = 23338001.40046683 #rad/s
Omegae = 0.0 #rad/s
Omegai = 0.0 #rad/s

k_nyquist = 3.140064641598747 #1/cm
kalias_em = 5.449663615983119 #1/cm
normk = 2.098480398345973 #1/cm
normfreq = 1000000000.0 #rad/s

kprime(k) = asin(k * dx/2.)*2./dx #Finite grid effects
oprime(o) = 2./dt * asin(dt/2. * o) #Finite grid effects
#kprime(k) = k
#oprime(o) = o

set yrange [0:8.]
set xrange [-0.15:0.15]
set cbrange [1e-2:1e5]
set samples nx

plot "run0_1d3v/disp_E2_x.dat" matrix using (($1-nx/2)*2*pi/dx/nx/normk/n_lineout_raumschritte):($2*2*pi/dt/normfreq/n_timesteps):($3*n_lineout_raumschritte*n_lineout_raumschritte) with image notitle \
, 2/(dt * normfreq) * asin(sqrt(dt**2/4.0 * (4*c**2 * (sin(normk*dx/2* x))**2 / dx**2 + omegape**2 + omegapi**2))) lt 3 lc rgb "black" title "EM mode" \
, sqrt(omegape**2 + omegapi**2)/normfreq lt 2 lc rgb "black" title "plasma freq {/Symbol w}_{p}" \
