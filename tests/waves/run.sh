#!/bin/bash
set -e


dir="run0_1d3v";
mkdir $dir;
time ../../python/1d3v/pypic.py --mime 1836.0 --Nt 10400 --Nx 2048 --nppc 256 --vthe 0.05 --outputdir $dir --rescale_dt 1.22638766024344700058 --rescale_dx 1.49895 --vbeamx 0.0 --vbeamy 0.0 --vbeamz 0.0 --deltav 0.0 --m 0 --theta 0.5 --Bx 0.0 --By 0.0 --Bz 0.0 --outputsteps 0 --particlesteps 0 | tee $dir/run.log
./run0_disp_E1_x.gp
./run0_disp_E2_x.gp

dir="run1_1d1v_noB";
mkdir $dir;
time ../../python/1d1v/pypic_noB.py --mime 1836.0 --Nt 10400 --Nx 2048 --nppc 128 --vthe 0.05 --outputdir $dir --rescale_dt 1.22638766024344700058 --rescale_dx 1.49895 --vbeam 0.0 --deltav 0.0 --m 0 --theta 0.5 --outputsteps 0 --particlesteps 0 | tee $dir/run.log
./run1_noB_disp_E0_x.gp


dir="run1_1d1v";
mkdir $dir;
time ../../python/1d1v/pypic.py --mime 1836.0 --Nt 10400 --Nx 2048 --nppc 128 --vthe 0.05 --outputdir $dir --rescale_dt 1.22638766024344700058 --rescale_dx 1.49895 --vbeam 0.0 --deltav 0.0 --m 0 --theta 0.5 --Bx 0. --outputsteps 0 --particlesteps 0 | tee $dir/run.log
./run1_disp_E0_x.gp


dir="run2_1d3v";
mkdir $dir;
time ../../python/1d3v/pypic.py --mime 1836.0 --Nt 10400 --Nx 2048 --nppc 256 --vthe 0.05 --outputdir $dir --rescale_dt 1.22638766024344700058 --rescale_dx 1.49895 --vbeamx 0.0 --vbeamy 0.0 --vbeamz 0.0 --deltav 0.0 --m 0 --theta 0.5 --Bx 28.43 --By 0.0 --Bz 0.0 --outputsteps 0 --particlesteps 0 | tee $dir/run.log
./run2_disp_E1_x.gp
./run2_disp_El_x.gp
./run2_disp_Er_x.gp


dir="run3_1d3v";
mkdir $dir;
time ../../python/1d3v/pypic.py --mime 1836.0 --Nt 10400 --Nx 2048 --nppc 256 --vthe 0.05 --outputdir $dir --rescale_dt 1.22638766024344700058 --rescale_dx 1.49895 --vbeamx 0.0 --vbeamy 0.0 --vbeamz 0.0 --deltav 0.0 --m 0 --theta 0.5 --Bx 0.0 --By 28.43 --Bz 0.0 --outputsteps 0 --particlesteps 0 --tol 1e-5 | tee $dir/run.log
./run3_disp_E0_x.gp
./run3_disp_E1_x.gp
./run3_disp_E2_x.gp


dir="run4_1d3v";
mkdir $dir;
time ../../python/1d3v/pypic.py --mime 1836.0 --Nt 2560 --Nx 512 --nppc 128 --vthe 0.05 --outputdir $dir --rescale_dt 1.0 --rescale_dx 0.093684375 --vbeamx 0.0 --vbeamy 0.0 --vbeamz 0.0 --deltav 0.0 --m 0 --theta 0.5 --Bx 28.43 --By 0.0 --Bz 0.0 --outputsteps 0 --particlesteps 0 --tol 1e-5 | tee $dir/run.log
./run4_disp_E1_x.gp
./run4_disp_El_x.gp
./run4_disp_Er_x.gp


dir="run5_1d3v";
mkdir $dir;
time ../../python/1d3v/pypic.py --mime 18.36 --Nt 3534 --Nx 3415 --Nspecies 2 --nppc 128 --vthe 0.05 --outputdir $dir --rescale_dt 0.1 --rescale_dx .3124 --vbeamx 0.0 --vbeamy 0.0 --vbeamz 0.0 --deltav 0.0 --m 0 --theta 0.5 --Bx 28.4282 --By 0. --Bz 0. --outputsteps 0 --particlesteps 0
./run5_disp_El_x.gp
./run5_disp_Er_x.gp
