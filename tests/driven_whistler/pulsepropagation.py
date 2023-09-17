#!/usr/bin/python3

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

c = 2.9979e10

# commandline
parser = argparse.ArgumentParser(description='Check charge conservation')
parser.add_argument('--indir', dest='indir', action='store', default="run0_1d3v", type=str, help="Input directory")
args = parser.parse_args()

f = h5py.File(args.indir+"/dings.h5", 'r')
# background magnetic field
B0  = np.array(f["/Timestep_0/felder/B/B[0]"])
Bg = np.mean(B0)

# find timesteps
times = []
for group in f:
    if group.startswith("Timestep_"):
        t = group.split("_")[1]
        times.append(int(t))

# transverse magnetic field
fout = open(args.indir+"/B2.dat", 'w')
fout.write("#B0 = "+str(Bg)+"\n")
for t in sorted(times):
    # remove ghost cells and unresolved dimensions
    B2  = np.array(f["/Timestep_"+str(t)+"/felder/B/B[2]"])[4:-4,0,0]
    fout.write(" ".join(map(str, B2))+"\n")
fout.close()
f.close()
