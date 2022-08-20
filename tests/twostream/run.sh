#!/bin/bash
set -e


mkdir 1d1v_noB
time ../../python/1d1v/pypic_noB.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir 1d1v_noB --rescale_dt .55530483123306872000 --vbeam 0.2 --deltav 0.001 --m 3 --rescale_dx .10185916357881301489 --theta 0.5 | tee 1d1v_noB/run.log

./phasespace.py --infile 1d1v_noB/electron/0.h5 --title r'$t = 0 \Delta{t}$' --outfile 1d1v_noB/phasespace_0.png
./phasespace.py --infile 1d1v_noB/electron/1000.h5 --title r'$t = 1000 \Delta{t}$' --outfile 1d1v_noB/phasespace_1000.png

mkdir 1d1v
time ../../python/1d1v/pypic.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir 1d1v --rescale_dt .55530483123306872000 --vbeam 0.2 --deltav 0.001 --m 3 --rescale_dx .10185916357881301489 --theta 0.5 --Bx 0. | tee 1d1v/run.log

./phasespace.py --infile 1d1v/electron/0.h5 --title r'$t = 0 \Delta{t}$' --outfile 1d1v/phasespace_0.png
./phasespace.py --infile 1d1v/electron/1000.h5 --title r'$t = 1000 \Delta{t}$' --outfile 1d1v/phasespace_1000.png


mkdir 1d3v
time ../../python/1d3v/pypic.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir 1d3v --rescale_dt .55530483123306872000 --vbeamx 0.2 --deltav 0.001 --m 3 --rescale_dx .10185916357881301489 --theta 0.5 --Bx 0. | tee 1d3v/run.log

./phasespace.py --infile 1d3v/electron/0.h5 --title r'$t = 0 \Delta{t}$' --outfile 1d3v/phasespace_0.png
./phasespace.py --infile 1d3v/electron/1000.h5 --title r'$t = 1000 \Delta{t}$' --outfile 1d3v/phasespace_1000.png



./plot_fig1.gp
./plot_fig2.gp
./plot_fig3.gp
