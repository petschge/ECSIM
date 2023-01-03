# Energy Conserving Semi-Implicit Particle-in-Cell (ECSIM)

This code is a test implementation of the semi-implicit energy conserving
Particle-in-Cell (PIC) algorithm proposed by G. Lapenta [1].

## Versions of the Code

Currently three versions of the code are implemented in python:

1. 1d1v phase space, only solving for one component of electric field, $E_{x}$.
   The magnetic field $B_{x}$ is a constant, since it is required to be
   homogeneous due to the divergence constraint.

2. 1d1v phase space, but retaining both $E_{x}$ and $B_{x}$. While $B_{x}$
   still needs to remain homogeneous and static, it is included in the Maxwell
   solver as a dynamical variable. This version is used to verify that the Maxwell
   solver (the most complicated part of the algorithm) preserves constant $B_x$.

3. 1d3v phase space, retaining all six components of the electromagnetic field.

All three version of the code make use of [Numpy](https://numpy.org/) and
[Numba](https://numba.pydata.org/) to improve the run time and use
[GMRES](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html)
from [Scipy](https://scipy.org/) to solve the linear system of equations that
gives the electromagnetic fields at the next timestep.

## Test Problems

The following test problems are currently implemented and are supplied with
simple scripts that run the code and produce plots:

### Two-stream instability

Two counter-streaming electron populations are set up. The diagnostic plots
include the growth rate of the electrostatic two-stream instability,
conservation of total energy, evolution of particle momentum and implied
effective numerical collision frequency.

### Weibel instability

Two counter-streaming electron populations are set up. Since the drift speed is
in a direction perpendicular to the resolved dimension, only the 1d3v code can
simulate this instability. The diagnostic plots include the growth rate of the
electromagnetic Weibel instability and conservation of total energy.

### Finite grid instability

Shows the absence of the finite grid instability even when increasingly
under-resolving the Debye length. Conservation of total energy and absence of
self-heating is shown. Fluctuations in kinetic and electric energy _decrease_
with larger cells and numerical collisionality drops.

### Plasma Wave Modes

This test follows the approach suggested in [2]. Eigenmodes in thermal
plasma under several different conditions are extracted by online Fourier
transform of the fields at run time and are compared to analytic
expectations from cold plasma theory.

### Thermalization test

This tests sets up hot protons and cold electrons ($T_{i}$ / $T_{e}$ = $10^{4}$,
$m_{i} / m_{e} = 100$) and looks at the heating rate of the electrons for
different particles per cell and different super-sampling factor $\Xi$. The
general trend is that moderate $\Xi$ in the range 2...100 heat electrons faster
than a simulation that resolves electron scales. At very large $\Xi$ the
heating rate drops significantly. Comparison with the explicit particle-in-cell
code VPIC [3] shows that VPIC heats the electrons faster at the same number of
particles per cell, which is expected due to the lower order of interpolation.
At a larger number of particles per cell which required comparable
computational effort due to the faster particles pushes in VPIC, the heating
rate in VPIC is lower thanin ECSIM.

### Acknowledgments

The code was developed during a technology evaluation phase for a project sponsored by US National Science Foundation under award NSF-2031024.

###References

[1] G. Lapenta, Exactly energy conserving semi-implicit particle in cell formulation, Journal of Computational Physics 334 (2017), p.349–366, [https://dx.doi.org/10.1016/j.jcp.2017.01.002](https://dx.doi.org/10.1016/j.jcp.2017.01.002)

[2] Patrick Kilian, Patricio A. Muñoz, Cedric Schreiner, and Felix Spanier, Plasma waves as a benchmark problem, Journal of Plasma Physics, 83(1), 707830101, [https://doi.org/10.1017/S0022377817000149](https://doi.org/10.1017/S0022377817000149)

[3] K.J. Bowers, B.J. Albright, B. Bergen and T.J.T. Kwan, Ultrahigh performance three-dimensional electromagnetic relativistic kinetic plasma simulation, Phys. Plasmas 15, 055703 (2008); [https://dx.doi.org/10.1063/1.2840133](https://dx.doi.org/10.1063/1.2840133)
