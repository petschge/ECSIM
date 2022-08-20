#!/usr/bin/python3

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

c = 2.9979e10

# commandline
parser = argparse.ArgumentParser(description='Check charge conservation')
parser.add_argument('--indir', dest='indir', action='store', default="1d1v_Xi1", type=str, help="Input directory")
args = parser.parse_args()

# Xi value
parts = args.indir.split("_")
if parts[-1].startswith("Xi"):
    Xi = float(parts[-1][2:])
else:
    Xi = 1.

# final thermal speed
f = h5py.File(args.indir+"/electron/1000.h5", 'r')
vx  = np.array(f["vx"])
f.close()
vth = np.std(vx)

# total energy
data = np.loadtxt(args.indir+"/energy.dat")
Etotal = data[:,1]
Etotal_err = np.average((Etotal-Etotal[0])/ Etotal[0])
if Etotal_err == 0.:
    Etotal_err = 1e-17

# electric energy
EE = data[:,5]
EE_err = np.std(EE/Etotal[0])
if EE_err == 0.:
    EE_err = 1e-17

# particle energy
Ekin = data[:,2]
Ekin_err = np.std(Ekin/Etotal[0])
if Ekin_err == 0.:
    Ekin_err = 1e-17

# collisionality
data2 = np.loadtxt(args.indir+"/momentum.dat")
t  = data2[:,0]
px = data2[:,2]
p0 = data2[0,1]
nu = (px[1:]-px[:-1]) / (t[1:]-t[:-1]) / p0
#plt.plot(t,px)
#plt.plot(0.5*(t[1:]+t[:1]),nu)
#plt.show()
nu_average=np.std(nu)

# output
print(Xi,vth/c,Etotal_err,EE_err,Ekin_err,nu_average)
