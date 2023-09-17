#!/usr/bin/python3

import numpy as np
import sys
from scipy.sparse.linalg import gmres
from timeit import default_timer as timer
import h5py
import argparse
from numba import njit
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter("error")

def ensure_dir(dirname):
    """
    Ensure that a directory with the supplied name exists
    """
    import os
    try:
        os.makedirs(dirname)
    except FileExistsError:
        pass


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Energy-Conserving Implicit PiC in 1d3v')
parser.add_argument('--mime', dest='mime', action='store', default=18.362, type=float, help="Mass ratio m_i/m_e (default: %(default)s)")
parser.add_argument('--vthe', dest='vthe', action='store', default=0.05, type=float, help="v_{th,e}/c (default: %(default)s)")
parser.add_argument('--TiTe', dest='TiTe', action='store', default=1., type=float, help="T_i / T_e (default: %(default)s)")
parser.add_argument('--wpe', dest='wpe', action='store', default=1e9, type=float, help="w_pe (default: %(default)s)")
parser.add_argument('--vbeamx', dest='vbeamx', action='store', default=0.0, type=float, help="v_{beam,x}/c (default: %(default)s)")
parser.add_argument('--vbeamy', dest='vbeamy', action='store', default=0.0, type=float, help="v_{beam,y}/c (default: %(default)s)")
parser.add_argument('--vbeamz', dest='vbeamz', action='store', default=0.0, type=float, help="v_{beam,z}/c (default: %(default)s)")
parser.add_argument('--deltav', dest='deltav', action='store', default=0.0, type=float, help="deltav/c (default: %(default)s)")
parser.add_argument('--m', dest='m', action='store', default=3, type=int, help="unstable mode m (default: %(default)s)")
parser.add_argument('--Nt', dest='Nt', action='store', default=1, type=int, help="number of time steps (default: %(default)s)")
parser.add_argument('--Nx', dest='Nxinner', action='store', default=2048, type=int, help="number of cells (default: %(default)s)")
parser.add_argument('--Nspecies', dest='Nspecies', action='store', default=1, type=int, help="number of species (default: %(default)s)")
parser.add_argument('--presolve', dest='presolve', action='store_true', help="Perform Poisson solve for E_x before timestepping (default: %(default)s)")
parser.add_argument('--theta', dest='theta', action='store', default=0.5, type=float, help="theta (default: %(default)s)")
parser.add_argument('--tol', dest='tol', action='store', default=1e-6, type=float, help="tol (default: %(default)s)")
parser.add_argument('--atol', dest='atol', action='store', default=1e-15, type=float, help="atol (default: %(default)s)")
parser.add_argument('--Bx', dest='B0x', action='store', default=28.43, type=float, help="Bx (default: %(default)s)")
parser.add_argument('--By', dest='B0y', action='store', default=0., type=float, help="By (default: %(default)s)")
parser.add_argument('--Bz', dest='B0z', action='store', default=0., type=float, help="Bz (default: %(default)s)")
parser.add_argument('--nppc', dest='nppc', action='store', default=64, type=int, help="number of particles per cell and species (default: %(default)s)")
parser.add_argument('--outputdir', dest='outputdir', action='store', default="run1", type=str, help="Output directory (default: %(default)s)")
parser.add_argument('--rescale_dx', dest='rescale_dx', action='store', default=1.0, type=float, help="Spatial resolution improvement (default: %(default)s)")
parser.add_argument('--rescale_dt', dest='rescale_dt', action='store', default=1.0, type=float, help="Temporal resolution improvement (default: %(default)s)")
parser.add_argument('--outputsteps', dest='outputsteps', action='store', default=1, type=int, help="Timesteps between field outputs (default: %(default)s)")
parser.add_argument('--particlesteps', dest='particlesteps', action='store', default=1000, type=int, help="Timesteps between particle outputs (default: %(default)s)")
parser.add_argument('--notiming', dest='timing', action='store_false', help="Suppress timing info")
parser.add_argument('--antenna_a0', dest='antenna_a0', action='store', default=0., type=float, help="strength of the driven current")
parser.add_argument('--antenna_w0', dest='antenna_w0', action='store', default=0., type=float, help="Driving frequency of antenna")
parser.add_argument('--antenna_delta_t', dest='antenna_delta_t', action='store', default=0., type=float, help="Duration of wave packet")
parser.add_argument('--antenna_t0', dest='antenna_t0', action='store', default=0., type=float, help="when the envelope is max")
parser.add_argument('--antenna_dw_dt', dest='antenna_dw_dt', action='store', default=0., type=float, help="frequency sweep rate")
parser.add_argument('--antenna_delta_x', dest='antenna_delta_x', action='store', default=0., type=float, help="source region")

args = parser.parse_args()

class my_timers:
    def __init__(self):
        self.timers = {
            'init':0.,
            'misc1':0.,
            'misc2':0.,
            'misc3':0.,
            'misc4':0.,
            'misc5':0.,
            'misc6':0.,
            'alpha':0.,
            'build_A':0.,
            'build_b':0.,
            'current':0.,
            'density':0.,
            'energies':0.,
            'gmres':0.,
            'localB':0.,
            'localE':0.,
            'mass_matrix':0.,
            'maxwell':0.,
            'newv':0.,
            'newx':0.,
            'timestepping':0.,
            'total':0.,
            }
        return

    def tic(self, _timer_):
        self.timers[_timer_] -= timer()
        return

    def toc(self, _timer_):
        self.timers[_timer_] += timer()
        return


Ng = 4 # number of ghost cells

args.Nxouter = args.Nxinner + 2*Ng
Np = args.Nxinner * args.nppc
np.random.seed(42)



c  = 2.9979e10    # cm/s
e  = 4.8032e-10   # statcoul
me = 9.1094e-28   # g
kB = 1.3807e-16   # erg/K
kBeV = 1.6022e-12 # erg/eV






timers = my_timers()
timers.tic("total")
timers.tic("init")

class Species:
    """
    One species of particles, such as electrons or protons or ions. Could be a subspecies such as hot electrons or a beam population.
    """

    def __init__(self, name, q, m, macro, Np, vth):
        """
        Create particles at t=0

        Np (int): number of particles to create
        vth (float): Isotropic thermal speed in cm/s

        Creates: ndarray(Np,float) for x in cm
           and 3 ndarray(Np,float) for vx,vy,vz in cm/s

        Return: nothing
        """

        self.name = name
        self.q = q
        self.m = m
        self.macro = macro

        #self.x = np.random.uniform(low=(Ng*dx), high=((args.Nxinner+Ng)*dx), size=Np)
        self.x = np.zeros(Np)
        nppc = Np//args.Nxinner
        p = 0
        for i in range(Ng,args.Nxinner+Ng):
            for j in range(nppc):
                self.x[p] = (i+np.random.uniform(0.,1.))*dx
                p += 1

        if args.vbeamx == 0. and args.vbeamy == 0. and args.vbeamz == 0.:
            # thermal plasma
            self.vx = np.random.normal(scale=vth, size=Np)
            self.vy = np.random.normal(scale=vth, size=Np)
            self.vz = np.random.normal(scale=vth, size=Np)
        elif args.vbeamx != 0. and args.vbeamy == 0. and args.vbeamz == 0.:
            # Two-stream instability
            self.vx = np.zeros(Np)
            Lx = args.Nxinner*dx
            for i in range(Np):
                self.vx[i] = (-1.)**i * args.vbeamx*c + np.random.normal(scale=vth, size=1) + args.deltav*c * np.sin(2.*np.pi*args.m*self.x[i]/Lx)
            self.vy = np.random.normal(scale=vth, size=Np)
            self.vz = np.random.normal(scale=vth, size=Np)
        elif args.vbeamx == 0. and args.vbeamy != 0. and args.vbeamz == 0.:
            # Weibel instability
            self.vx = np.random.normal(scale=vth, size=Np)
            self.vy = np.zeros(Np)
            Lx = args.Nxinner*dx
            for i in range(Np):
                self.vy[i] = (-1.)**i * args.vbeamy*c + np.random.normal(scale=vth, size=1) + args.deltav*c * np.sin(2.*np.pi*args.m*self.x[i]/Lx)
            self.vz = np.random.normal(scale=vth, size=Np)
        else:
            sys.stderr.write("This combination of beam velocities does not correspond to a known test problem.\n")
            sys.stderr.write("You will have to provide your own particle setup in "+sys.argv[0]+"\n")
            sys.exit(1)


    def uncenter(self):
        """
        Pushes particle positions back from t=0 to t^-1/2

        Reads: x ndarray(Np,float): position in cm
               vx,vy,vz ndarray(Np,float): velocities in cm/s

        Changes: x to ndarray(Np,float) for oldx in cm

        Returns: nothing
        """
        Np = len(self.x)
        assert len(self.vx) == Np
        assert len(self.vy) == Np
        assert len(self.vz) == Np

        # push x back from t^1/2 to t=0
        oldx = self.x - self.vx * dt/2.

        # periodic wrap around
        cond = oldx < Ng*dx
        oldx[cond] += args.Nxinner*dx

        cond = oldx > (args.Nxinner+Ng)*dx
        oldx[cond] -= args.Nxinner*dx

        self.x = oldx

    def newx(self):
        """
        Update position x from t^n-1/2 to t^n+1/2 using v at t^n

        Reads: x ndarray(Np,float): position in cm
               vx,vy,vz ndarray(Np,float): velocity in cm/s

        Changes: x to ndarray(Np,float) for new positions as floats in cm, periodically
                 wrapped around if the particle leaves the (interior of) the domain

        Returns: nothing
        """

        timers.tic("newx")
        newx = self.x + self.vx*dt

        cond = newx < Ng*dx
        newx[cond] += args.Nxinner*dx

        cond = newx > (args.Nxinner+Ng)*dx
        newx[cond] -= args.Nxinner*dx

        self.x = newx
        timers.toc("newx")

    def newv(self, Elocalx,Elocaly,Elocalz):
        """
        Update velocity v from t^n to t^n+1 using E^n+theta at position x at t^n+1/2

        Elocalx,Elocaly,Elocalz ndarray(Np,float): electric field in statV/cm

        Reads: x ndarray(Np,float): position in cm
               vx,vy,vz ndarray(Np,float): velocity in cm/s
               alpha ndarray((3,3,Np),float): dimensionless matrix for B field rotation

        Changes: ndarray(Np,float) for new positions as floats in cm, periodically
                 wrapped around if the particle leaves the (interior of) the domain

        Return: nothing
        """

        timers.tic("newv")

        Np = len(self.x)
        assert len(self.vx) == Np
        assert len(self.vy) == Np
        assert len(self.vz) == Np
        assert self.alpha.shape[2] == Np
        assert len(Elocalx) == Np
        assert len(Elocaly) == Np
        assert len(Elocalz) == Np

        beta = (self.q*dt)/(2.*self.m)

        vrx = self.vx + beta*Elocalx
        vry = self.vy + beta*Elocaly
        vrz = self.vz + beta*Elocalz

        # alpha is dimensionless
        vpx = self.alpha[0,0]*vrx + self.alpha[0,1]*vry + self.alpha[0,2]*vrz
        vpy = self.alpha[1,0]*vrx + self.alpha[1,1]*vry + self.alpha[1,2]*vrz
        vpz = self.alpha[2,0]*vrx + self.alpha[2,1]*vry + self.alpha[2,2]*vrz

        assert len(vpx) == Np
        assert len(vpy) == Np
        assert len(vpz) == Np

        self.vx = 2.*vpx - self.vx
        self.vy = 2.*vpy - self.vy
        self.vz = 2.*vpz - self.vz

        timers.toc("newv")

    def get_localB(self, Bx,By,Bz):
        """
        Interpolate magnetic field from face centers to position of each particle

        Bx,By,Bz ndarray(Nxouter,float): Magnetic field in Gauss

        Reads: x ndarray(Np,float): Particle position in cm

        Returns: 3 ndarray(Np,float) for the local magnetic field in Gauss at each
                 particle location
        """
        timers.tic("localB")

        assert len(Bx) == args.Nxouter
        assert len(By) == args.Nxouter
        assert len(Bz) == args.Nxouter

        cellx = np.array((self.x/dx), dtype=int)

        wxhl = W_CIC(self.x, cellx-1, 0.5)
        wxhc = W_CIC(self.x, cellx  , 0.5)
        wxhr = W_CIC(self.x, cellx+1, 0.5)

        Blocalx = wxhl * Bx[cellx-1] + wxhc * Bx[cellx  ] + wxhr * Bx[cellx+1]
        Blocaly = wxhl * By[cellx-1] + wxhc * By[cellx  ] + wxhr * By[cellx+1]
        Blocalz = wxhl * Bz[cellx-1] + wxhc * Bz[cellx  ] + wxhr * Bz[cellx+1]

        timers.toc("localB")
        return Blocalx,Blocaly,Blocalz

    def get_localE(self, Ex,Ey,Ez):
        """
        Interpolate electric field from nodes to position of each particle

        Ex,Ey,Ez ndarray(Nxouter,float): Electric field in statV/cm

        Reads: x ndarray(Np,float): Particle position in cm

        Returns: 3 ndarray(Np,float) for the local electric field in statV/cm at
                 each particle location
        """
        timers.tic("localE")

        assert len(Ex) == args.Nxouter
        assert len(Ey) == args.Nxouter
        assert len(Ez) == args.Nxouter

        cellx = np.array((self.x/dx), dtype=int)

        wxfl = W_CIC(self.x, cellx-1, 0.0)
        wxfc = W_CIC(self.x, cellx  , 0.0)
        wxfr = W_CIC(self.x, cellx+1, 0.0)

        Elocalx = wxfl * Ex[cellx-1] + wxfc * Ex[cellx  ] + wxfr * Ex[cellx+1]
        Elocaly = wxfl * Ey[cellx-1] + wxfc * Ey[cellx  ] + wxfr * Ey[cellx+1]
        Elocalz = wxfl * Ez[cellx-1] + wxfc * Ez[cellx  ] + wxfr * Ez[cellx+1]

        timers.toc("localE")
        return Elocalx,Elocaly,Elocalz

@njit
def compute_alpha(x, Bx,By,Bz, q,m):
    """
    Compute alpha matrix representing rotation due to the magnetic field

    x ndarray(Np,float): position in cm
    Bx,By,Bz ndarray(Np,float): magnetic field at particle locations in Gauss
    q float: charge of the particle in statC
    m float: mass of the particle in g

    Changes: alpha ndarray((3,3,Np), float) for the dimensionless matrix for each particle
    """

    Np = len(x)
    assert len(Bx) == Np
    assert len(By) == Np
    assert len(Bz) == Np

    beta = (q*dt)/(2.*m*c) # inverse units of B
    bx = (beta*Bx) # dimensionless
    by = (beta*By)
    bz = (beta*Bz)

    # Papers since at least Vu, H. X., & Brackbill, J. U. (1992) define this as
    # (\mathbb{1} - \beta\mathbb{1}\times B + \beta^2 B B) / (1+\beta^2 B^2)
    # abusing the notation of the cross product to be between a matrix and a vector
    # and writing the outer product of two vectors without operator symbol

    # if this was transposed, it would basically flips the sign of the charge
    # This way round is correct, because the transpose gives left hand circular whistlers
    alpha = np.zeros( (3,3,Np) )
    alpha[0,0,:] = ( 1. + bx*bx) / (1. + bx*bx + by*by + bz*bz)
    alpha[0,1,:] = ( bz + bx*by) / (1. + bx*bx + by*by + bz*bz)
    alpha[0,2,:] = (-by + bx*bz) / (1. + bx*bx + by*by + bz*bz)
    alpha[1,0,:] = (-bz + bx*by) / (1. + bx*bx + by*by + bz*bz)
    alpha[1,1,:] = ( 1. + by*by) / (1. + bx*bx + by*by + bz*bz)
    alpha[1,2,:] = ( bx + by*bz) / (1. + bx*bx + by*by + bz*bz)
    alpha[2,0,:] = ( by + bx*bz) / (1. + bx*bx + by*by + bz*bz)
    alpha[2,1,:] = (-bx + by*bz) / (1. + bx*bx + by*by + bz*bz)
    alpha[2,2,:] = ( 1. + bz*bz) / (1. + bx*bx + by*by + bz*bz)

    #alpha = np.array( [ [ 1. + bx*bx,  bz + bx*by, -by + bx*bz],
    #                    [-bz + bx*by,  1. + by*by,  bx + by*bz],
    #                    [ by + bx*bz, -bx + by*bz,  1. + bz*bz] ] ) / (1. + bx*bx + by*by + bz*bz)

    # between 1. and 1 - 1.4e-8. but systematically below 1.
    # # alpha should have determiante 1 since it is just a rotation. let's check
    # det = np.zeros(Np)
    # for p in range(Np):
    #     det[p] = np.linalg.det(alpha[:,:,p])
    #     #alpha[:,:,p] /= np.cbrt(det[p]) # this makes it worse
    # amin = np.amin(det)
    # amax = np.amax(det)
    # print("Entries in alpha are between "+str(amin)+" and "+str(amax))

    return alpha

@njit
def density(Nx, x, q, macro):
    """
    Charge deposition

    Nx int: Nxouter
    x ndarray(Np,float): particle position in cm at t^n+1/2

    Returns: ndarray(Nxouter,float) for charge at t^n+1/2 in statC/cm^3 already
             reduced to the interior of the periodic domain at cell centers
    """
    rho = np.zeros(Nx)

    q = q * macro / dx**3

    for p in range(len(x)):
        cellx = int(x[p]/dx)

        for i in range(-2,3):
            w = W_CIC_scalar_numba(x[p], cellx+i, 0.5)
            rho[cellx+i] += q * w

    rho = charge_boundary(rho)
    return rho

@njit
def ecsim_current(Nx, x, vx,vy,vz, alpha, q,macro):
    """
    Current deposition of \hat{J}_{sg}

    Nx int: Nxouter
    x ndarray(Np,float): particle position in cm at t^n+1/2
    vx,vy,vz ndarray(Np,float): particle velocities in cm/s at t^n
    alpha ndarray((2,2,Np),float): dimensionless rotation matrixes

    Returns: 3 ndarray(Nxouter,float) for currents at t^n in statC/cm^2/s already
             reduced to the interior of the periodic domain
    FIXME: is that current really at t^n+1/2?
    FIXME: I assume j is at integer options same as E?
    """
    jx,jy,jz = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)

    q = q * macro / dx**3

    Np = len(x)
    assert len(vx) == Np
    assert len(vy) == Np
    assert len(vz) == Np
    assert alpha.shape[2] == Np

    for p in range(Np):
        cellx = int(x[p]/dx)

        #for i in range(-2,3):
        for i in range(-1,2): # this seems fine
            w = W_CIC_scalar_numba(x[p], cellx+i, 0.0)
            jx[cellx+i] += (q * w * (alpha[0,0,p]*vx[p] + alpha[0,1,p]*vy[p] + alpha[0,2,p]*vz[p]))
            jy[cellx+i] += (q * w * (alpha[1,0,p]*vx[p] + alpha[1,1,p]*vy[p] + alpha[1,2,p]*vz[p]))
            jz[cellx+i] += (q * w * (alpha[2,0,p]*vx[p] + alpha[2,1,p]*vy[p] + alpha[2,2,p]*vz[p]))


    jx = charge_boundary(jx)
    jy = charge_boundary(jy)
    jz = charge_boundary(jz)
    return jx,jy,jz

@njit
def W_CIC_scalar_numba(xp, xg, deltayee):
    """
    Integral of the zeroth order tophat-shaped cloud-in-cell shape function

    xp (float): Particle position in cm
    xg (int): Grid index in cells
    deltayee(float): Either 0. if the grid index is supposed to be at an
                     integer node position
                     or 0.5 if the grid index referes to a quantity at
                     half-integer position at a face center

    Returns: weight w as a dimensionless float
    """

    x = np.abs((xp/dx) - xg - deltayee)
    if x<1.0:
        w = 1.-x
    else:
        w = 0.

    return w

@njit
def matrix_boundary(M):
    """
    Move contribution that ended in any of the ghost zones off a matrix to the
    periodically wrapped interior of the domain

    M: ndarray((Nxouter,Nxouter), float)

    Returns: numpy array where ghost zones are zeroed out and contributions
             moved to the interior
    """
    M[-2*Ng: -Ng,:] += M[   :Ng,:]
    M[   Ng:2*Ng,:] += M[-Ng:  ,:]
    M[     :Ng  ,:]  = 0. * M[0,0]
    M[  -Ng:    ,:]  = 0. * M[0,0]

    M[:,-2*Ng: -Ng] += M[:,   :Ng]
    M[:,   Ng:2*Ng] += M[:,-Ng:  ]
    M[:,     :Ng  ]  = 0. * M[0,0]
    M[:,  -Ng:    ]  = 0. * M[0,0]

    return M

@njit
def naive_numba_mass_matrix(Nx, x, alpha, q, m, macro, M00,M01,M02, M10,M11,M12, M20,M21,M22):
    """
    Compute mass matrices that contain the reaction of the magnetised plasma to
    the electric field

    Nx int: Nxouter
    x ndarray(Np,float): particle positions in cm at t^n+1/2
    alpha ndarray((3,3,Np), float): dimensionless rotation matrices at each
                                    particle location to describe B field

    Returns: 9 ndarray((Nxouter,Nxouter),float) for the mass matrices. The units of each
             entry are 1/seconds since we have to multiply an electric field in
             statV/cm and get a current in statC/cm^2/s
    """
    #assert len(x) == Np
    #assert alpha.shape[2] == Np
    assert M00.shape[0] == Nx
    assert M00.shape[1] == Nx
    assert M01.shape[0] == Nx
    assert M01.shape[1] == Nx
    assert M02.shape[0] == Nx
    assert M02.shape[1] == Nx
    assert M10.shape[0] == Nx
    assert M10.shape[1] == Nx
    assert M11.shape[0] == Nx
    assert M11.shape[1] == Nx
    assert M12.shape[0] == Nx
    assert M12.shape[1] == Nx
    assert M20.shape[0] == Nx
    assert M20.shape[1] == Nx
    assert M21.shape[0] == Nx
    assert M21.shape[1] == Nx
    assert M22.shape[0] == Nx
    assert M22.shape[1] == Nx

    Vg = dx**3
    beta = (q*dt)/(2.*m*c) # dimension l^1/2 t m^-1/2

    for p in range(np.int64(len(x))):
        cellx = int(x[p]/dx)
        for i in range(-2,3):
            for j in range(-2,3):
                w = beta/Vg * q*macro*c * W_CIC_scalar_numba(x[p], cellx+i, 0.0) * W_CIC_scalar_numba(x[p], cellx+j, 0.0)
                M00[cellx+i,cellx+j] += w * alpha[0,0,p]
                M01[cellx+i,cellx+j] += w * alpha[0,1,p]
                M02[cellx+i,cellx+j] += w * alpha[0,2,p]
                M10[cellx+i,cellx+j] += w * alpha[1,0,p]
                M11[cellx+i,cellx+j] += w * alpha[1,1,p]
                M12[cellx+i,cellx+j] += w * alpha[1,2,p]
                M20[cellx+i,cellx+j] += w * alpha[2,0,p]
                M21[cellx+i,cellx+j] += w * alpha[2,1,p]
                M22[cellx+i,cellx+j] += w * alpha[2,2,p]

    M00 = matrix_boundary(M00)
    M01 = matrix_boundary(M01)
    M02 = matrix_boundary(M02)
    M10 = matrix_boundary(M10)
    M11 = matrix_boundary(M11)
    M12 = matrix_boundary(M12)
    M20 = matrix_boundary(M20)
    M21 = matrix_boundary(M21)
    M22 = matrix_boundary(M22)

    #return M00,M01,M02,M10,M11,M12,M20,M21,M22

def init_fields(B0x,B0y,B0z):
    """
    Create fields at t=0

    B0x, B0y, B0z (float): Homogeneous background magnetic field in Gauss

    Returns: 3 ndarray(Nxouter,float) for Bx,By,Bz in Gauss
         and 3 ndarray(Nxouter,float) for Ex,Ey,Ez in statC/cm
    """

    Bx = np.full(  args.Nxouter, B0x)
    By = np.full(  args.Nxouter, B0y)
    Bz = np.full(  args.Nxouter, B0z)
    Ex = np.zeros( args.Nxouter )
    Ey = np.zeros( args.Nxouter )
    Ez = np.zeros( args.Nxouter )
    return Bx,By,Bz, Ex,Ey,Ez

def field_boundary(V):
    """
    Fills ghost zones on the left and right with values from the right/left
    interior of the domain

    V: ndarray(Nxouter,float)

    Returns: numpy array where 0:Ng and Nxinner+Ng:Nxouter are filled in
    """
    V[   :Ng] = V[-2*Ng: -Ng]
    V[-Ng:  ] = V[   Ng:2*Ng]
    return V

@njit
def charge_boundary(rho):
    """
    Move contribution that ended in the ghost zones on the left and right to
    the interior of the domain on the right/left and zero out the ghost zones.

    rho: ndarray(Nxouter,float)

    Returns: numpy array where 0:Ng and Nxinner+Ng:Nxouter are zeroed out and contributions moved to the interior
    """
    rho[-2*Ng: -Ng] += rho[   :Ng]
    rho[   Ng:2*Ng] += rho[-Ng:  ]
    rho[     :Ng  ]  = 0. * rho[0]
    rho[  -Ng:    ]  = 0. * rho[0]

    return rho

def W_CIC(xp, xg, deltayee):
    """
    Integral of the zeroth order tophat-shaped cloud-in-cell shape function

    xp (float): Particle position in cm
    xg (int): Grid index in cells
    deltayee(float): Either 0. if the grid index is supposed to be at an
                     integer node position
                     or 0.5 if the grid index referes to a quantity at
                     half-integer position at a face center

    Calling with arrays ndarray(Np,float), ndarray(Np,int), ndarray(Np,float) works too

    Returns: weight w as a dimensionless float
    """

    w = np.zeros_like(xp, dtype=float)
    x = np.abs((xp/dx) - xg - deltayee)
    cond = x<1.0
    w[cond] = 1.0 - x[cond]

    return w

def W_TSC(xp, xg, deltayee):
    """
    Integral of the first order triangular-shaped cloud shape function

    xp (float): Particle position in cm
    xg (int): Grid index in cells
    deltayee(float): Either 0. if the grid index is supposed to be at an
                     integer node position
                     or 0.5 if the grid index referes to a quantity at
                     half-integer position at a face center

    Calling with arrays ndarray(Np,float), ndarray(Np,int), ndarray(Np,float) works too

    Returns: weight w as a dimensionless float
    """

    w = np.zeros_like(xp, dtype=float)
    x = np.abs(xp/dx - xg - deltayee)
    cond = x<0.5
    w[cond] = 0.75 - x[cond]*x[cond]
    cond2 = np.logical_and(np.logical_not(cond), x<1.5)
    w[cond2] = 0.5 * (1.5-x[cond2])**2

    return w

def poisson(rho):
    rhok = np.fft.fft(rho[Ng:-Ng])
    kx = np.fft.fftfreq(args.Nxinner, d=dx/(2.*np.pi))
    for i in range(rhok.shape[0]):
        if kx[i] == 0.:
            rhok[i] = 0. # assume neutralizing background
        else:
            rhok[i] *= 4. * np.pi / (kx[i]**2)
    phi = np.fft.ifft(rhok).real
    phiout = np.zeros( args.Nxouter )
    phiout[Ng:args.Nxinner+Ng] = phi
    phiout = field_boundary(phiout)
    return phiout

def ES_field(phi):
    Ex = np.zeros(args.Nxouter)

    for i in range(1, args.Nxouter-1):
        Ex[i] = - (phi[i] - phi[i-1]) / dx

    Ex = field_boundary(Ex)
    return Ex

def maxwell(Ex,Ey,Ez, Bx,By,Bz, jhatx,jhaty,jhatz, M00,M01,M02,M10,M11,M12,M20,M21,M22, theta):
    """
    Solve linear problem to get updated electric and magnetic field

    Ex,Ey,Ez ndarray(Nxouter,float): electric field at t^n in statV/cm
    Bx,By,Bz ndarray(Nxouter,float): magnetic field at t^n in Gauss

    jhatx,jhaty,jhatz ndarray(Nxouter,float): explicit current in statC/cm^2/s

    M00,M01,M02,M10,M11,M12,M20,M21,M22 ndarray((Nxouter,Nxouter),float): mass matricies in 1/s

    Returns: 3 ndarray(Nxouter,float) for E at t^n+1 in statV/cm
             3 ndarray(Nxouter,float) for B at t^n+1 in Gauss
    """
    timers.tic("maxwell")

    assert len(Ex) == args.Nxouter
    assert len(Ey) == args.Nxouter
    assert len(Ez) == args.Nxouter
    assert len(Bx) == args.Nxouter
    assert len(By) == args.Nxouter
    assert len(Bz) == args.Nxouter
    assert len(jhatx) == args.Nxouter
    assert len(jhaty) == args.Nxouter
    assert len(jhatz) == args.Nxouter
    assert M00.shape[0] == args.Nxouter
    assert M00.shape[1] == args.Nxouter
    assert M01.shape[0] == args.Nxouter
    assert M01.shape[1] == args.Nxouter
    assert M02.shape[0] == args.Nxouter
    assert M02.shape[1] == args.Nxouter
    assert M10.shape[0] == args.Nxouter
    assert M10.shape[1] == args.Nxouter
    assert M11.shape[0] == args.Nxouter
    assert M11.shape[1] == args.Nxouter
    assert M12.shape[0] == args.Nxouter
    assert M12.shape[1] == args.Nxouter
    assert M20.shape[0] == args.Nxouter
    assert M20.shape[1] == args.Nxouter
    assert M21.shape[0] == args.Nxouter
    assert M21.shape[1] == args.Nxouter
    assert M22.shape[0] == args.Nxouter
    assert M22.shape[1] == args.Nxouter

    # form combined vector for fields. they have compatible units of m^1/2 l^-1/2 t^-1 in CGS
    F = np.concatenate( (Ex[Ng:-Ng],Ey[Ng:-Ng],Ez[Ng:-Ng],Bx[Ng:-Ng],By[Ng:-Ng],Bz[Ng:-Ng]) )

    timers.tic("build_b")
    # our right hand side has units of Gauss/(cm/s * s) = (statV/cm)/cm = m^1/2 l^-3/2 t^-1
    b = np.zeros(6*(args.Nxinner))

    jx = jhatx + (M00@Ex + M01@Ey + M02@Ez)*(1.-theta)
    jy = jhaty + (M10@Ex + M11@Ey + M12@Ez)*(1.-theta)
    jz = jhatz + (M20@Ex + M21@Ey + M22@Ez)*(1.-theta)

    for m in range(0*args.Nxinner,1*args.Nxinner):
        i = m - 0*args.Nxinner
        b[m] = (- F[0*args.Nxinner+i] / (c*dt) + 4*np.pi/c * jx[i+Ng])
    for m in range(1*args.Nxinner,2*args.Nxinner):
        i = m - 1*args.Nxinner
        iprev = i-1
        if iprev < 0:
            iprev += args.Nxinner
        b[m] = ((theta-1.)/dx * (-F[5*args.Nxinner+i] + F[5*args.Nxinner+iprev]) - F[1*args.Nxinner+i] / (c*dt) + 4*np.pi/c * jy[i+Ng])
    for m in range(2*args.Nxinner,3*args.Nxinner):
        i = m - 2*args.Nxinner
        iprev = i-1
        if iprev < 0:
            iprev += args.Nxinner
        b[m] = ((theta-1.)/dx * (F[4*args.Nxinner+i] - F[4*args.Nxinner+iprev]) - F[2*args.Nxinner+i] / (c*dt) + 4*np.pi/c * jz[i+Ng])

    for m in range(3*args.Nxinner,4*args.Nxinner):
        i = m - 3*args.Nxinner
        b[m] = (F[3*args.Nxinner+i] / (c*dt))
    for m in range(4*args.Nxinner,5*args.Nxinner):
        i = m - 4*args.Nxinner
        inext = i+1
        if inext >= args.Nxinner:
            inext -= args.Nxinner
        b[m] = ((theta-1.)/dx * (-F[2*args.Nxinner+inext] + F[2*args.Nxinner+i]) + F[4*args.Nxinner+i] / (c*dt))
    for m in range(5*args.Nxinner,6*args.Nxinner):
        i = m - 5*args.Nxinner
        inext = i+1
        if inext >= args.Nxinner:
            inext -= args.Nxinner
        b[m] = ((theta-1.)/dx * (F[1*args.Nxinner+inext] - F[1*args.Nxinner+i]) + F[5*args.Nxinner+i] / (c*dt))

    timers.toc("build_b")

    # are are constructing a matrix M for the equation
    # M.F = b where
    # F has units of      m^1/2 l^-1/2 t^-1
    # and b has units of  m^1/2 l^-3/2 t^-1
    # this implies that all entries in M should have units 1/l
    timers.tic("build_A")
    A = np.zeros( (6*args.Nxinner,6*args.Nxinner) )

    A[0*args.Nxinner:1*args.Nxinner, 0*args.Nxinner:1*args.Nxinner] -= (M00[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    A[0*args.Nxinner:1*args.Nxinner, 1*args.Nxinner:2*args.Nxinner] -= (M01[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    A[0*args.Nxinner:1*args.Nxinner, 2*args.Nxinner:3*args.Nxinner] -= (M02[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    for m in range(0*args.Nxinner,1*args.Nxinner):
        i = m - 0*args.Nxinner
        A[m,i] += (-1./(c*dt))

    A[1*args.Nxinner:2*args.Nxinner, 0*args.Nxinner:1*args.Nxinner] -= (M10[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    A[1*args.Nxinner:2*args.Nxinner, 1*args.Nxinner:2*args.Nxinner] -= (M11[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    A[1*args.Nxinner:2*args.Nxinner, 2*args.Nxinner:3*args.Nxinner] -= (M12[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    for m in range(1*args.Nxinner,2*args.Nxinner):
        i = m - 1*args.Nxinner
        iprev = i-1
        if iprev < 0:
            iprev += args.Nxinner
        A[m,5*args.Nxinner+i    ] += (-theta/dx  )
        A[m,5*args.Nxinner+iprev] += ( theta/dx  )
        A[m,  args.Nxinner+i    ] += ( -1./(c*dt))

    A[2*args.Nxinner:3*args.Nxinner, 0*args.Nxinner:1*args.Nxinner] -= (M20[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    A[2*args.Nxinner:3*args.Nxinner, 1*args.Nxinner:2*args.Nxinner] -= (M21[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    A[2*args.Nxinner:3*args.Nxinner, 2*args.Nxinner:3*args.Nxinner] -= (M22[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    for m in range(2*args.Nxinner,3*args.Nxinner):
        i = m - 2*args.Nxinner
        iprev = i-1
        if iprev < 0:
            iprev += args.Nxinner
        A[m,4*args.Nxinner+i    ] += ( theta/dx )
        A[m,4*args.Nxinner+iprev] += (-theta/dx )
        A[m,2*args.Nxinner+i    ] += (-1./(c*dt))

    for m in range(3*args.Nxinner,4*args.Nxinner):
        i = m - 3*args.Nxinner
        A[m,3*args.Nxinner+i] += (1./(c*dt))
    for m in range(4*args.Nxinner,5*args.Nxinner):
        i = m - 4*args.Nxinner
        inext = i+1
        if inext >= args.Nxinner:
            inext -= args.Nxinner
        A[m,2*args.Nxinner+inext] += (-theta/dx)
        A[m,2*args.Nxinner+i    ] += ( theta/dx)
        A[m,4*args.Nxinner+i    ] += (1./(c*dt))
    for m in range(5*args.Nxinner,6*args.Nxinner):
        i = m - 5*args.Nxinner
        inext = i+1
        if inext >= args.Nxinner:
            inext -= args.Nxinner
        A[m,  args.Nxinner+inext] += ( theta/dx)
        A[m,  args.Nxinner+i    ] += (-theta/dx)
        A[m,5*args.Nxinner+i    ] += (1./(c*dt))

    timers.toc("build_A")

    timers.tic("gmres")
    info = 1
    tol = args.tol
    atol = args.atol
    while info > 0 and tol < 1e-2 and atol < 1e-2:
        Fnext,info = gmres(A,b, tol=tol/c, atol=atol, x0=F)
        if info > 0:
            tol *= 10
            atol *= 10
    timers.toc("gmres")

    if info < 0:
        sys.stderr.write("Illegal input in timestep "+str(t)+"\n")
        sys.exit(1)
    elif info > 0:
        sys.stderr.write("Did not convergence in timestep "+str(t)+"\n")
        sys.exit(1)

    newEx = np.zeros(args.Nxouter)
    newEy = np.zeros(args.Nxouter)
    newEz = np.zeros(args.Nxouter)
    newBx = np.zeros(args.Nxouter)
    newBy = np.zeros(args.Nxouter)
    newBz = np.zeros(args.Nxouter)
    newEx[Ng:-Ng] = Fnext[0*args.Nxinner:1*args.Nxinner]
    newEy[Ng:-Ng] = Fnext[1*args.Nxinner:2*args.Nxinner]
    newEz[Ng:-Ng] = Fnext[2*args.Nxinner:3*args.Nxinner]
    newBx[Ng:-Ng] = Fnext[3*args.Nxinner:4*args.Nxinner]
    newBy[Ng:-Ng] = Fnext[4*args.Nxinner:5*args.Nxinner]
    newBz[Ng:-Ng] = Fnext[5*args.Nxinner:6*args.Nxinner]

    newEx = field_boundary(newEx)
    newEy = field_boundary(newEy)
    newEz = field_boundary(newEz)
    newBx = field_boundary(newBx)
    newBy = field_boundary(newBy)
    newBz = field_boundary(newBz)

    timers.toc("maxwell")

    return newEx,newEy,newEz, newBx,newBy,newBz

def saveh5(A, name):
    f = h5py.File(args.outputdir+"/dings.h5", 'a')
    parts = name.split("/")
    grp = f
    for p in parts[1:-1]:
        try:
            grp = grp.create_group(p)
            B = np.array([dx])
            grp.create_dataset("dx", data=B)
            B = np.array([dt])
            grp.create_dataset("dt", data=B)
        except ValueError:
            grp = grp[p]
    if len(A.shape) == 0:
        shape = (1,1,1)
    elif len(A.shape) == 1:
        shape = (A.shape[0],1,1)
    elif len(A.shape) == 2:
        shape = (A.shape[0],A.shape[1],1)
    elif len(A.shape) == 3:
        shape = (A.shape[0],A.shape[1],A.shape[2])
    else:
        sys.stderr.write("Can save array with "+str(len(A.shape))+" dimensions to "+name+"\n")
        sys.exit(1)
    # FIXME: should we write units as an attribute?
    dset = grp.create_dataset(parts[-1], data=A, shape=shape)
    f.close()

def save1d(A, name, scale):
    f = open(args.outputdir+"/"+name, 'w')
    if scale is None:
        scale = np.max(np.abs(A))
    if scale <= 0.:
        scale = 1.
    for i in range(0, A.shape[0]):
        f.write(str(i)+" "+str(A[i]/scale)+"\n")
    f.close()

def save2d(A, name, scale):
    f = open(args.outputdir+"/"+name, 'w')
    if scale is None:
        scale = np.max(np.abs(A))
    if scale <= 0.:
        scale = 1.
    for j in range(0, A.shape[1]):
        f.write(" ".join(map(str,A[:,j]/scale)) + "\n")
    f.close()

def save_particles(s,t):
    dirname = args.outputdir+"/"+s.name
    ensure_dir(dirname)

    h5f = h5py.File(dirname+"/"+str(t)+".h5", 'w')
    h5f.create_dataset("x", data=s.x)
    h5f.create_dataset("vx", data=s.vx)
    h5f.create_dataset("vy", data=s.vy)
    h5f.create_dataset("vz", data=s.vz)
    h5f.close()

def save_energies(energies, t, Bx,By,Bz, Ex,Ey,Ez, species):
    """
    Save enery at time t^n+1

    energies: open file pointer
    t: int time step we just completed
    Ex: electric field at t^n+1
    species: list of species that each has: position t^n+1/2 and velocity t^n+1
    """
    timers.tic("energies")

    eBx = (np.sum(Bx[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eBy = (np.sum(By[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eBz = (np.sum(Bz[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eB  = eBx + eBy + eBz
    eEx = (np.sum(Ex[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eEy = (np.sum(Ey[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eEz = (np.sum(Ez[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eE  = eEx + eEy + eEz
    eFields = eB + eE
    eParticles = np.zeros(len(species))
    for i,s in enumerate(species):
        E = 0.5 * s.m * s.macro * np.sum(s.vx**2 + s.vy**2 + s.vz**2)
        eParticles[i] = E

    energies.write(str((t*dt*args.wpe)))                   # 1
    energies.write(" "+str((np.sum(eParticles)+eFields)))  # 2
    energies.write(" "+str(np.sum(eParticles)))            # 3
    energies.write(" "+str(eFields))                       # 4
    energies.write(" "+str(eB))                            # 5
    energies.write(" "+str(eE))                            # 6
    energies.write(" "+str(eBx))                           # 7
    energies.write(" "+str(eBy))
    energies.write(" "+str(eBz))                           # 9
    energies.write(" "+str(eEx))
    energies.write(" "+str(eEy))                           # 11
    energies.write(" "+str(eEz))                           # 12
    for i in range(len(eParticles)):
        energies.write(" "+str(eParticles[i]))
    energies.write("\n")
    energies.flush()


    timers.toc("energies")



#
# main
#
ensure_dir(args.outputdir)

vthe = args.vthe * c
sys.stderr.write("wpe = "+str(args.wpe)+"rad/s\n")
ne = args.wpe**2 * me / (4.*np.pi * e**2)
sys.stderr.write("n0 = "+str(ne)+"cm^{-3}\n")

sys.stderr.write("vthe = "+str(vthe)+"cm/s\n")
sys.stderr.write("     = "+str((vthe/c))+"c\n")
Te = vthe**2*me/kB
sys.stderr.write("Te = "+str(Te)+"K\n")
TeeV = vthe**2*me/kBeV
sys.stderr.write("   = "+str(TeeV)+"eV\n")

lD = vthe/args.wpe
de = c/args.wpe
sys.stderr.write("lD = "+str(lD)+"cm\n")
sys.stderr.write("   = "+str(lD/de)+"de\n")
sys.stderr.write("de = "+str(de)+"cm\n")

dx = lD/args.rescale_dx
sys.stderr.write("dx = "+str(dx)+"cm\n")
sys.stderr.write("   = "+str(dx/lD)+"l_D\n")
sys.stderr.write("   = "+str(dx/de)+"d_e\n")

Lx = dx*args.Nxinner
sys.stderr.write("Lx = "+str(Lx)+"cm\n")
sys.stderr.write("   = "+str(Lx/lD)+"l_D\n")
sys.stderr.write("   = "+str(Lx/de)+"d_e\n")

B0 = np.sqrt(args.B0x**2 + args.B0y**2 + args.B0z**2)
Wce = e*B0/(me*c)
sys.stderr.write("B0 = "+str(args.B0x)+", "+str(args.B0y)+","+str(args.B0z)+" G\n")
sys.stderr.write("|B0| = "+str(B0)+"G\n")

dt_wpe = 1./args.wpe
dt_Wce = 1./Wce if Wce>0. else np.inf
dt_part = dx/(4.*vthe)
dt_yee = 0.9999 * dx/(np.sqrt(2.)*c)
dt = min(dt_wpe,dt_Wce,dt_part,dt_yee) / args.rescale_dt

sys.stderr.write("dt = "+str(dt)+"s\n")
sys.stderr.write("wpe*dt = "+str((args.wpe*dt))+"\n")
sys.stderr.write("Wce*dt = "+str((Wce*dt))+"\n")
sys.stderr.write("c*dt/dx = "+str((c*dt/dx))+"\n")
sys.stderr.write("wpe/Wce = "+str(args.wpe/Wce if Wce>0. else np.inf)+"\n")
sys.stderr.write("beta_e = "+str(8.*np.pi*kB*Te/B0**2 if B0>0. else np.inf)+"\n")

sys.stderr.write("T = "+str(dt*args.Nt)+"s\n")
sys.stderr.write("  = "+str(dt*args.Nt*args.wpe)+"wpe^-1\n")
sys.stderr.write("  = "+str(dt*args.Nt*Wce)+"Wce^-1\n")

if args.Nspecies > 1:
    wpi = args.wpe / np.sqrt(args.mime)
    Wci = Wce / args.mime
    di = c/wpi
    sys.stderr.write("  = "+str(dt*args.Nt*wpi)+"wpi^-1\n")
    sys.stderr.write("  = "+str(dt*args.Nt*Wci)+"Wci^-1\n")
    sys.stderr.write("wpi*dt = "+str((wpi*dt))+"\n")
    sys.stderr.write("Wci*dt = "+str((Wci*dt))+"\n")
    sys.stderr.write("di = "+str(di)+"cm\n")
    sys.stderr.write("dx = "+str(dx/di)+"d_i\n")
    sys.stderr.write("Lx = "+str(Lx/di)+"d_i\n")

    vA = B0/np.sqrt(4.*np.pi*ne*me*args.mime) if B0 > 0. else 0.
    sys.stderr.write("vA = "+str(vA)+" cm/s\n")
    sys.stderr.write("   = "+str(vA/c)+" c\n")
    if B0 > 0.:
        sys.stderr.write("Wci/vA = "+str(Wci/vA)+" 1/cm\n")
        sys.stderr.write("vA/Wci = "+str(vA/Wci)+" cm\n")


sys.stderr.write("Nt = "+str(args.Nt)+"\n")

Nphys = ne * dx**3
macro = Nphys / args.nppc
sys.stderr.write("macro = "+str(macro)+"\n")

if args.antenna_a0 > 0. and args.antenna_w0 > 0. and args.antenna_delta_t > 0. and args.antenna_t0 > 0. and args.antenna_delta_x > 0.:
    # kmin resolved in the domain
    kmin = 2.*np.pi / Lx
    # desired k0
    k0_prime = np.sqrt(args.antenna_w0/Wce) / de
    # actual k0 such that it's an even number of kmin
    args.antenna_k0 = int(0.5*k0_prime/kmin)*2*kmin
    # actual w0 that matches that k0
    args.antenna_w0 = (args.antenna_k0*de)**2 * Wce

    sys.stderr.write("a0 = "+str(args.antenna_a0)+"not sure of the units\n")
    sys.stderr.write("w0 = "+str(args.antenna_w0)+"rad/s\n")
    sys.stderr.write("   = "+str(args.antenna_w0/args.wpe)+"wpe\n")
    sys.stderr.write("   = "+str(args.antenna_w0/Wce)+"Wce\n")
    sys.stderr.write("k0 = "+str(args.antenna_k0)+"m^-1\n")
    sys.stderr.write("   = "+str(args.antenna_k0*dx)+"/dx\n")
    sys.stderr.write("   = "+str(args.antenna_k0*de)+"/de\n")
    sys.stderr.write("   = "+str(args.antenna_k0*Lx)+"/Lx\n")
    sys.stderr.write("   = "+str(args.antenna_k0*dx/(2.*np.pi))+"*2pi/dx\n")
    sys.stderr.write("   = "+str(args.antenna_k0*de/(2.*np.pi))+"*2pi/de\n")
    sys.stderr.write("   = "+str(args.antenna_k0*Lx/(2.*np.pi))+"*2pi/Lx\n")
    sys.stderr.write("delta_t = "+str(args.antenna_delta_t)+"Hz\n")
    sys.stderr.write("        = "+str(args.antenna_delta_t*Wce)+"Wce^-1\n")
    sys.stderr.write("        = "+str(args.antenna_delta_t*args.antenna_w0)+"w0^-1\n")
    sys.stderr.write("        = "+str(args.antenna_delta_t*args.wpe)+"wpe^-1\n")
    sys.stderr.write("t0 = "+str(args.antenna_t0)+"Hz\n")
    sys.stderr.write("   = "+str(args.antenna_t0/args.antenna_delta_t)+"delta_t\n")
    sys.stderr.write("   = "+str(args.antenna_t0*Wce)+"Wce^-1\n")
    sys.stderr.write("   = "+str(args.antenna_t0*args.antenna_w0)+"w0^-1\n")
    sys.stderr.write("   = "+str(args.antenna_t0*args.wpe)+"wpe^-1\n")
    sys.stderr.write("dw/dt = "+str(args.antenna_dw_dt)+"rad/s^2\n")
    sys.stderr.write("      = "+str(args.antenna_dw_dt*args.antenna_delta_t/args.antenna_w0)+"w0/delta_t\n")
    sys.stderr.write("delta_x = "+str(args.antenna_delta_x)+"m\n")
    sys.stderr.write("        = "+str(args.antenna_delta_x*args.antenna_k0)+"/k0\n")


t=0

Bx,By,Bz, Ex,Ey,Ez = init_fields(args.B0x,args.B0y,args.B0z)

species = []
if args.Nspecies > 0:
    electrons = Species("electron", -e, me,           macro, Np, vthe)
    species.append(electrons)
if args.Nspecies > 1:
    protons   = Species("proton",    e, args.mime*me, macro, Np, vthe*np.sqrt(args.TiTe/args.mime))
    species.append(protons)
if args.Nspecies > 2:
    sys.stderr.write("Please define what other species you want\n")
    sys.exit(1)

# get density at t^0
rho = np.zeros(args.Nxouter)
for s in species:
    timers.tic('density')
    rhos = density(args.Nxouter, s.x, s.q, s.macro)
    rho += rhos
    timers.toc('density')

Q0 = np.sum(rho)

if args.presolve:
    # solve for Ex at t^0
    phi = poisson(rho)
    Ex = ES_field(phi)

# compute error charge rho-(divE)/4pi
rhoE = rho - (np.roll(Ex,-1)-Ex)/dx / (4.*np.pi)
rhoE = field_boundary(rhoE)

for s in species:
    # Pushes particle positions back from t=0 to t^-1/2
    s.uncenter()

# f = open(args.outputdir+"/traj.dat", "w")

energies = open(args.outputdir+"/energy.dat", "w")
dispEx = np.zeros( (args.Nt, args.Nxinner) )
dispEy = np.zeros( (args.Nt, args.Nxinner) )
dispEz = np.zeros( (args.Nt, args.Nxinner) )

# f.write("0 "+str(electrons.x[0])+" "+str(electrons.vx[0])+" "+str(electrons.vy[0])+" "+str(electrons.vz[0])+"\n")

# reset HDF5 output
h5f = h5py.File(args.outputdir+"/dings.h5", 'w')
cfg_grp = h5f.create_group("config")
cfg_grp.create_dataset("plasmafreq", data=np.array([args.wpe]))
cfg_grp.create_dataset("B0", data=np.array([args.B0x,args.B0y,args.B0z]))
cfg_grp.create_dataset("units", data=np.array("cgs", dtype=h5py.string_dtype('ascii', 3)))
catfile = "lineout_raumschritte = 1\nmpzume="+str(args.mime)+"\nwidthbg="+str(vthe/c)+"\n"
utf8_type=h5py.string_dtype('utf-8', len(catfile))
catbytes = np.array(catfile.encode("utf-8"), dtype=utf8_type)
cfg_grp.create_dataset("cat-file", data=catbytes)
h5f.close()

saveh5(rho,  "/Timestep_0/rho/rhoL/rhoL[0]")
saveh5(rhoE, "/Timestep_0/rho/rhoE/rhoE")
saveh5(Ex, "/Timestep_0/felder/E/E[0]")
saveh5(Ey, "/Timestep_0/felder/E/E[1]")
saveh5(Ez, "/Timestep_0/felder/E/E[2]")
saveh5(Bx, "/Timestep_0/felder/B/B[0]")
saveh5(By, "/Timestep_0/felder/B/B[1]")
saveh5(Bz, "/Timestep_0/felder/B/B[2]")
if args.particlesteps > 0:
    for s in species:
        save_particles(s, 0)

energies.write("#t ")
energies.write(" total")
energies.write(" particles")
energies.write(" fields")
energies.write(" B")
energies.write(" E")
energies.write(" Bx")
energies.write(" By")
energies.write(" Bz")
energies.write(" Ex")
energies.write(" Ey")
energies.write(" Ez")
for s in species:
    energies.write(" "+s.name)
energies.write("\n")
energies.write("# dt*wpe = "+str(dt*args.wpe)+"\n")
energies.write("# dt*Wce = "+str(dt*Wce)+"\n")
energies.flush()
save_energies(energies, 0, Bx,By,Bz, Ex,Ey,Ez, species)

energycheck = open(args.outputdir+"/energycheck.dat", 'w')
charge = open(args.outputdir+"/charge.dat", 'w')

momentum = open(args.outputdir+"/momentum.dat", "w")
p0 = 0.
px,py,pz = np.zeros(args.Nt+1),np.zeros(args.Nt+1),np.zeros(args.Nt+1)
for s in species:
    p0 += s.m * s.macro * args.vbeamx*c * len(s.vx)
    px[0] += s.m * s.macro * np.sum(s.vx)
    py[0] += s.m * s.macro * np.sum(s.vy)
    pz[0] += s.m * s.macro * np.sum(s.vz)
if args.vbeamx == 0.:
    for s in species:
        p0 += me * s.macro * args.vthe*c * len(s.vx)
momentum.write(str(0.)+" "+str(p0)+" "+str(px[0])+" "+str(py[0])+" "+str(pz[0])+"\n")

timers.toc("init")
timers.tic("timestepping")

for t in range(1,args.Nt+1):
    timers.tic("misc1")
    # net current from all species
    jx,jy,jz = np.zeros(args.Nxouter), np.zeros(args.Nxouter), np.zeros(args.Nxouter)

    # total mass matrix for all species
    M00,M01,M02 = np.zeros( (args.Nxouter,args.Nxouter) ), np.zeros( (args.Nxouter,args.Nxouter) ), np.zeros( (args.Nxouter,args.Nxouter) )
    M10,M11,M12 = np.zeros( (args.Nxouter,args.Nxouter) ), np.zeros( (args.Nxouter,args.Nxouter) ), np.zeros( (args.Nxouter,args.Nxouter) )
    M20,M21,M22 = np.zeros( (args.Nxouter,args.Nxouter) ), np.zeros( (args.Nxouter,args.Nxouter) ), np.zeros( (args.Nxouter,args.Nxouter) )
    timers.toc("misc1")

    for s in species:
        # Update x from t^n-1/2 to t^n+1/2 using v at t^n
        s.newx()

        # local magnetic fields
        Blocalx,Blocaly,Blocalz = s.get_localB(Bx,By,Bz)

        # compute alpha for all particles
        timers.tic("alpha")
        s.alpha = compute_alpha(s.x, Blocalx,Blocaly,Blocalz, s.q,s.m)
        timers.toc("alpha")

        # compute current \hat{j}_sg using x at t^n+1/2 and B and v at t^n
        timers.tic('current')
        jhatx,jhaty,jhatz =  ecsim_current(args.Nxouter, s.x, s.vx,s.vy,s.vz, s.alpha, s.q,s.macro)
        timers.toc('current')

        # save for analysis
        # if args.outputsteps > 0 and t%args.outputsteps == 0:
        #     saveh5(jhatx, "/Timestep_"+str(t)+"/"+s.name+"/jhat[1]")
        #     saveh5(jhaty, "/Timestep_"+str(t)+"/"+s.name+"/jhat[2]")
        #     saveh5(jhatz, "/Timestep_"+str(t)+"/"+s.name+"/jhat[3]")

        # add to total current
        timers.tic("misc2")
        jx += jhatx
        jy += jhaty
        jz += jhatz
        timers.toc("misc2")

        # compute mass matrices
        timers.tic('mass_matrix')
        naive_numba_mass_matrix(args.Nxouter, s.x, s.alpha, s.q, s.m, s.macro, M00,M01,M02,M10,M11,M12,M20,M21,M22)
        timers.toc('mass_matrix')

    # add antenna current
    if args.antenna_a0 > 0. and args.antenna_w0 > 0. and args.antenna_delta_t > 0. and args.antenna_t0 > 0. and args.antenna_delta_x > 0.:
        At = args.antenna_a0 * np.exp(-(t*dt-args.antenna_t0)**2 / args.antenna_delta_t**2)
        w = args.antenna_w0 + args.antenna_dw_dt*(t*dt - args.antenna_t0)
        for i in range(Ng, Ng+args.Nxinner):
            x = i*dx
            Ax = np.exp(-(x-0.5*Lx)**2 / args.antenna_delta_x**2)
            jy[i] += At*Ax * np.sin(args.antenna_k0*x - w*t*dt)
            jz[i] += At*Ax * np.cos(args.antenna_k0*x - w*t*dt)

    # save total current
    if args.outputsteps > 0 and t%args.outputsteps == 0:
        saveh5(jx,  "/Timestep_"+str(t)+"/rho/rhoL/rhoL[1]")
        saveh5(jy,  "/Timestep_"+str(t)+"/rho/rhoL/rhoL[2]")
        saveh5(jz,  "/Timestep_"+str(t)+"/rho/rhoL/rhoL[3]")

    # update E and B from t^n to t^n+1 using j at t^n+1/2 and the implicit theta algorithm
    newEx,newEy,newEz,newBx,newBy,newBz = maxwell(Ex,Ey,Ez, Bx,By,Bz, jx,jy,jz, M00,M01,M02,M10,M11,M12,M20,M21,M22, args.theta)

    # get electric field at t^n+1/2
    timers.tic("misc4")
    midEx = args.theta*newEx + (1.-args.theta)*Ex
    midEy = args.theta*newEy + (1.-args.theta)*Ey
    midEz = args.theta*newEz + (1.-args.theta)*Ez
    timers.toc("misc4")

    # save new fields
    if args.outputsteps > 0 and t%args.outputsteps == 0:
        saveh5(newEx, "/Timestep_"+str(t)+"/felder/E/E[0]")
        saveh5(newEy, "/Timestep_"+str(t)+"/felder/E/E[1]")
        saveh5(newEz, "/Timestep_"+str(t)+"/felder/E/E[2]")
        saveh5(newBx, "/Timestep_"+str(t)+"/felder/B/B[0]")
        saveh5(newBy, "/Timestep_"+str(t)+"/felder/B/B[1]")
        saveh5(newBz, "/Timestep_"+str(t)+"/felder/B/B[2]")

    jxgbar = np.zeros(args.Nxouter)
    jygbar = np.zeros(args.Nxouter)
    jzgbar = np.zeros(args.Nxouter)
    Ekinprev = 0.
    Ekinnext = 0.

    for s in species:
        # local electric field
        Elocalx,Elocaly,Elocalz = s.get_localE(midEx, midEy, midEz)

        # keep old velocity
        s.oldvx = s.vx.copy()
        s.oldvy = s.vy.copy()
        s.oldvz = s.vz.copy()

        # update particle velocities from t^n to t^n+1 using E at t^n+1
        s.newv(Elocalx,Elocaly,Elocalz)
        # check energy conservation
        s.vxbar = 0.5*(s.oldvx + s.vx)
        s.vybar = 0.5*(s.oldvy + s.vy)
        s.vzbar = 0.5*(s.oldvz + s.vz)
        jxsgbar,jysgbar,jzsgbar = ecsim_current(args.Nxouter, s.x, s.vxbar,s.vybar,s.vzbar, s.alpha, s.q,s.macro)
        jxgbar += jxsgbar
        jygbar += jysgbar
        jzgbar += jzsgbar
        Ekinprev += 0.5*s.m*s.macro*np.sum(s.oldvx**2)
        Ekinprev += 0.5*s.m*s.macro*np.sum(s.oldvy**2)
        Ekinprev += 0.5*s.m*s.macro*np.sum(s.oldvz**2)
        Ekinnext += 0.5*s.m*s.macro*np.sum(s.vx**2)
        Ekinnext += 0.5*s.m*s.macro*np.sum(s.vy**2)
        Ekinnext += 0.5*s.m*s.macro*np.sum(s.vz**2)

    deltaEparticle = Ekinnext-Ekinprev
    deltaEfieldsEx = (np.sum(newEx[Ng:-Ng]**2)/(8.*np.pi) - np.sum(Ex[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfieldsEy = (np.sum(newEy[Ng:-Ng]**2)/(8.*np.pi) - np.sum(Ey[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfieldsEz = (np.sum(newEz[Ng:-Ng]**2)/(8.*np.pi) - np.sum(Ez[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfieldsBx = (np.sum(newBx[Ng:-Ng]**2)/(8.*np.pi) - np.sum(Bx[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfieldsBy = (np.sum(newBy[Ng:-Ng]**2)/(8.*np.pi) - np.sum(By[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfieldsBz = (np.sum(newBz[Ng:-Ng]**2)/(8.*np.pi) - np.sum(Bz[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfields = deltaEfieldsEx+deltaEfieldsEy+deltaEfieldsEz + deltaEfieldsBx+deltaEfieldsBy+deltaEfieldsBz
    work = np.sum(jxgbar[Ng:-Ng]*midEx[Ng:-Ng] + jygbar[Ng:-Ng]*midEy[Ng:-Ng] + jzgbar[Ng:-Ng]*midEz[Ng:-Ng]) * dt * dx**3
    #print("Particle energy change: "+str(deltaEparticle)+", Work done: "+str(work)+", Field energy change:"+str(deltaEfields))
    energycheck.write(str(t)+" "+str(deltaEparticle)+" "+str(work)+" "+str(deltaEfields)+"\n")

    # get density at t^n+1/2
    rho = np.zeros(args.Nxouter)
    for s in species:
        timers.tic('density')
        rhos = density(args.Nxouter, s.x, s.q, s.macro)
        rho += rhos
        timers.toc('density')

    # check charge conservation
    Q = np.sum(rho)
    charge.write(str(t)+" "+str(Q/Q0)+"\n")

    # compute error charge rho-(divE)/4pi
    timers.tic("misc5")
    rhoE = rho - (np.roll(midEx,-1)-midEx)/dx / (4.*np.pi)
    rhoE = field_boundary(rhoE)
    timers.toc("misc5")

    if args.outputsteps > 0 and t%args.outputsteps == 0:
        saveh5(rho,  "/Timestep_"+str(t)+"/rho/rhoL/rhoL[0]")
        saveh5(rhoE, "/Timestep_"+str(t)+"/rho/rhoE/rhoE")

    # update electric and magnetic field to be at t^n+1
    timers.tic("misc6")
    Ex = newEx
    Ey = newEy
    Ez = newEz
    Bx = newBx
    By = newBy
    Bz = newBz
    timers.toc("misc6")

    if args.particlesteps > 0 and t%args.particlesteps == 0:
        for s in species:
            save_particles(s, t)

    # save energy at t^n+1
    save_energies(energies, t, Bx,By,Bz, Ex,Ey,Ez, species)

    dispEx[t-1,:] = Ex[Ng:-Ng] * np.sin(np.pi*t/args.Nt)**2
    dispEy[t-1,:] = Ey[Ng:-Ng] * np.sin(np.pi*t/args.Nt)**2
    dispEz[t-1,:] = Ez[Ng:-Ng] * np.sin(np.pi*t/args.Nt)**2

    # f.write(str(t)+" "+str(electrons.x[0])+" "+str(electrons.vx[0])+" "+str(electrons.vy[0])+" "+str(electrons.vz[0])+"\n")

    # plt.scatter(electrons.x/dx, electrons.vx/c)
    # plt.show()

    for s in species:
        px[t] += s.m * s.macro * np.sum(s.vx)
        py[t] += s.m * s.macro * np.sum(s.vy)
        pz[t] += s.m * s.macro * np.sum(s.vz)
    momentum.write(str(t*dt*args.wpe)+" "+str(p0)+" "+str(px[t])+" "+str(py[t])+" "+str(pz[t])+"\n")

disp_kw = np.fft.fft2(dispEx)
power = np.abs(np.fft.fftshift(disp_kw.T)[:,args.Nt//2:])
save2d(power, "disp_E0_x.dat", scale=-1.)
disp_kw = np.fft.fft2(dispEy)
power = np.abs(np.fft.fftshift(disp_kw.T)[:,args.Nt//2:])
save2d(power, "disp_E1_x.dat", scale=-1.)
disp_kw = np.fft.fft2(dispEz)
power = np.abs(np.fft.fftshift(disp_kw.T)[:,args.Nt//2:])
save2d(power, "disp_E2_x.dat", scale=-1.)

dispEl = dispEy - 1j * dispEz
disp_kw = np.fft.fft2(dispEl)
power = np.abs(np.fft.fftshift(disp_kw.T)[:,args.Nt//2:])
save2d(power, "disp_El_x.dat", scale=-1.)
dispEr = dispEy + 1j * dispEz
disp_kw = np.fft.fft2(dispEr)
power = np.abs(np.fft.fftshift(disp_kw.T)[:,args.Nt//2:])
save2d(power, "disp_Er_x.dat", scale=-1.)


# traj
# f.close()

timers.toc("timestepping")

energies.close()
momentum.close()

if p0 > 0.:
    nu = (px[1:] - px[:-1])/dt / p0
    nu_avg = np.std(nu/args.wpe)
    print("nu_avg = "+str(nu_avg)+" wpe")

timers.toc("total")

print("Total cost:     "+str(timers.timers["total"]))
if args.timing:
    print("Initalization:  "+str(timers.timers["init"]))
    print("Time Stepping:  "+str(timers.timers["timestepping"]))
    print("Mass Matrix:    "+str(timers.timers["mass_matrix"]))
    print("Maxwell cost:   "+str(timers.timers["maxwell"]))
    print("  build_b cost: "+str(timers.timers["build_b"]))
    print("  build_A cost: "+str(timers.timers["build_A"]))
    print("  GMRES cost:   "+str(timers.timers["gmres"]))
    print("Current cost:   "+str(timers.timers["current"]))
    print("Density cost:   "+str(timers.timers["density"]))
    print("Local B:        "+str(timers.timers["localB"]))
    print("Local E:        "+str(timers.timers["localE"]))
    print("New X:          "+str(timers.timers["newx"]))
    print("New V:          "+str(timers.timers["newv"]))
    print("Energies:       "+str(timers.timers["energies"]))
    print("Alpha:          "+str(timers.timers["alpha"]))
    print("Misc1:          "+str(timers.timers["misc1"]))
    print("Misc2:          "+str(timers.timers["misc2"]))
    print("Misc3:          "+str(timers.timers["misc3"]))
    print("Misc4:          "+str(timers.timers["misc4"]))
    print("Misc5:          "+str(timers.timers["misc5"]))
    print("Misc6:          "+str(timers.timers["misc6"]))

Np_total = 0
for s in species:
    Np_total += len(s.x)
print("particle updates per second: "+str((Np_total*args.Nt)/timers.timers["total"]))

energycheck.close()
charge.close()
