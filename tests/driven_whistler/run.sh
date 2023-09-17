#!/bin/bash
set -e

dir="run0_1d3v";
mkdir $dir;
time ../../python/1d3v/pypic.py                     \
    --mime 1836                                     \
    --vthe 0.02                                     \
    --wpe 56414.49647808635359452433                \
    --Nx 1000                                       \
    --rescale_dx 0.04                               \
    --Nspecies 2                                    \
    --Bx 0.000801875545471031                       \
    --antenna_w0 4231.08723585647625                \
    --antenna_delta_t .01181729357321531154         \
    --antenna_t0 .02                                \
    --antenna_dw_dt 35804.19839485519161503951      \
    --antenna_delta_x 10068557.89106258441322279363 \
    --antenna_a0 0.05                               \
    --Nt 10000                                      \
    --outputdir $dir

./pulsepropagation.py
./pulsepropagation.gp
./dispersionrelation.gp
./energy.gp
