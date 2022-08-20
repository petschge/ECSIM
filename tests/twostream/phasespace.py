#!/usr/bin/python3

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Check charge conservation')
parser.add_argument('--infile', dest='infile', action='store', default="run1", type=str, help="Input file")
parser.add_argument('--title', dest='title', action='store', default="", type=str, help="Plot title")
parser.add_argument('--outfile', dest='outfile', action='store', default="", type=str, help="filename for plot. interactive display if empty")
args = parser.parse_args()

f = h5py.File(args.infile, 'r')

x   = np.array(f["x"])
vx  = np.array(f["vx"])

f.close()

plt.scatter(x, vx/2.9979e10)
# plt.plot(X05, divE/(4.*np.pi), label="divE/4pi")
# plt.plot(X05, rho - divE/(4.*np.pi), label="rho - divE/4pi")
# plt.plot(X05, rhoE, label="rhoE")
# plt.legend()
plt.xlabel("x (cm)")
plt.ylabel("$v_x$ (c)")
plt.title(args.title)

if args.outfile == "":
    plt.show()
else:
    plt.savefig(args.outfile)
