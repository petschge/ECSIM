#!/bin/bash
set -e


mkdir 1d3v
time ../../python/1d3v/pypic.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir 1d3v --rescale_dt .55530483123306872000 --vbeamy 0.8 --deltav 0.001 --m 3 --rescale_dx .10185916357881301489 --theta 0.5 --Bx 0. --tol 1e-6 --atol 1e-14 --presolve | tee 1d3v/run.log

./phasespace.py --infile 1d3v/electron/0.h5 --title r'$t = 0 \Delta{t}$' --outfile 1d3v/phasespace_0.png
./phasespace.py --infile 1d3v/electron/1000.h5 --title r'$t = 1000 \Delta{t}$' --outfile 1d3v/phasespace_1000.png

./FFTmode.py
./plot_fig2.gp
