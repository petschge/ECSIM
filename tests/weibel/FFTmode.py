#!/usr/bin/python3

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Growth rate of the first few Fourier modes')
parser.add_argument('--inputdir', dest='inputdir', action='store', default="1d3v", type=str, help="Input directory")
args = parser.parse_args()

f = h5py.File(args.inputdir+"/dings.h5", 'r')
timesteps = []
for d in f.keys():
    if d.startswith("Timestep_"):
        parts = d.split("_")
        if len(parts) == 2:
            timesteps.append(int(parts[1]))
timesteps=sorted(timesteps)
Nt = len(timesteps)

times  = np.zeros(Nt)
energy1 = np.zeros(Nt)
energy2 = np.zeros(Nt)
energy3 = np.zeros(Nt)
energy4 = np.zeros(Nt)

wpe = f["/config/plasmafreq"]

for i,t in enumerate(timesteps):
    d = f["/Timestep_"+str(t)]
    dt = d["dt"][0]

    Bz   = np.array(d["felder"]["B"]["B[2]"])[:,0,0]
    kBz  = np.fft.fft(Bz)
    times[i]  = t*dt*wpe
    energy1[i] = abs(kBz[1])**2
    energy2[i] = abs(kBz[2])**2
    energy3[i] = abs(kBz[3])**2
    energy4[i] = abs(kBz[4])**2
    #energy[i] = np.sum(Bz**2)

f.close()

plt.plot(times, energy1, label="m=1", color='blue')
plt.plot(times, energy2, label="m=2", color='orange')
plt.plot(times, energy3, label="m=3", color='darkgray')
plt.plot(times, energy4, label="m=4", color='darkslategray')
gamma1 = 0.8 * np.sqrt((1**2)/(1+1**2))
plt.plot(times, 0.1*np.exp(times*gamma1*2), color='blue', linestyle="--", label="theory")
gamma2 = 0.8 * np.sqrt((2**2)/(1+2**2))
plt.plot(times, 0.1*np.exp(times*gamma2*2), color='orange', linestyle="--")
gamma3 = 0.8 * np.sqrt((3**2)/(1+3**2))
plt.plot(times, 0.1*np.exp(times*gamma3*2), color='darkgray', linestyle="--")
gamma4 = 0.8 * np.sqrt((4**2)/(1+4**2))
plt.plot(times, 0.1*np.exp(times*gamma4*2), color='darkslategray', linestyle="--")
plt.xlabel(r"$t \omega_{pe}$")
plt.xlim(0,50)
#plt.ylim(1e-5,1e6)
plt.ylim(1,1e6)
plt.legend()
plt.yscale('log')
#plt.show()
#plt.savefig(args.inputdir+"/figure1.png")
plt.savefig("figure1.png")
