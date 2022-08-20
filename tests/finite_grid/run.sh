#!/bin/bash
set -e


rm -f Xi_1d1v_noB.dat;
for logXi in $(seq 0 0.5 16); do
    Xi=$(echo "e(l(10)*$logXi) + 0.5" | bc -l | cut -d . -f 1);
    echo $Xi;
    dir="1d1v_noB_Xi"$Xi
    rescale_dx=$(echo ".10185916357881301489*e(l(10)*$logXi)" | bc -l);
    rescale_dt=$(echo ".55530483123306872000*e(l(10)*$logXi)" | bc -l);

    mkdir $dir;
    time ../../python/1d1v/pypic_noB.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir $dir --rescale_dt $rescale_dt --rescale_dx $rescale_dx --theta 0.5 | tee $dir/run.log
    ./analyze.py --indir $dir >> Xi_1d1v_noB.dat;
done


rm -f Xi_1d1v.dat;
for logXi in $(seq 0 0.5 16); do
    Xi=$(echo "e(l(10)*$logXi) + 0.5" | bc -l | cut -d . -f 1);
    echo $Xi;
    dir="1d1v_Xi"$Xi
    rescale_dx=$(echo ".10185916357881301489*e(l(10)*$logXi)" | bc -l);
    rescale_dt=$(echo ".55530483123306872000*e(l(10)*$logXi)" | bc -l);

    mkdir $dir;
    time ../../python/1d1v/pypic.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir $dir --rescale_dt $rescale_dt --rescale_dx $rescale_dx --theta 0.5 --Bx 0. | tee $dir/run.log
    ./analyze.py --indir $dir >> Xi_1d1v.dat;
done


rm -f Xi_1d3v.dat;
for logXi in $(seq 0 0.5 16); do
    Xi=$(echo "e(l(10)*$logXi) + 0.5" | bc -l | cut -d . -f 1);
    echo $Xi;
    dir="1d3v_Xi"$Xi
    rescale_dx=$(echo ".10185916357881301489*e(l(10)*$logXi)" | bc -l);
    rescale_dt=$(echo ".55530483123306872000*e(l(10)*$logXi)" | bc -l);

    mkdir $dir;
    time timeout 600 ../../python/1d3v/pypic.py --mime 1836.2 --Nt 1000 --Nx 64 --nppc 154 --vthe 0.01 --outputdir $dir --rescale_dt $rescale_dt --rescale_dx $rescale_dx --theta 0.5 --Bx 0. --tol 1e-5 --atol 1e-25 | tee $dir/run.log
    ./analyze.py --indir $dir >> Xi_1d3v.dat;
done


./plot_fig1.gp
