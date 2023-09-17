#!/bin/bash

for gen in 1 2 3; do
  nppc=$(echo "5*10^$gen" | bc -l)
  for xi in 1 3 5 10 30 100 300 1000; do
    ../../python/1d3v/pypic.py              \
      --mime 100                            \
      --vthe 0.005                          \
      --TiTe 1e4                            \
      --Nt $[1000000/$xi]                   \
      --Nx 32                               \
      --Nspecies 2                          \
      --theta 0.5                           \
      --Bx 0.003                            \
      --By 0.                               \
      --Bz 0.                               \
      --nppc ${nppc}                        \
      --outputdir thermal_equib${gen}_${xi} \
      --rescale_dx $(echo "1/$xi" | bc -l ) \
      --outputsteps 0                       \
      --particlesteps 0                     \
      --wpe 100000.;
  done;
done

./plot_energy.gp
./plot_Te.gp
