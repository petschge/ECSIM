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
parser.add_argument('--vbeam', dest='vbeam', action='store', default=0.0, type=float, help="v_{beam}/c (default: %(default)s)")
parser.add_argument('--deltav', dest='deltav', action='store', default=0.0, type=float, help="deltav/c (default: %(default)s)")
parser.add_argument('--m', dest='m', action='store', default=3, type=int, help="unstable mode m (default: %(default)s)")
parser.add_argument('--Nt', dest='Nt', action='store', default=1, type=int, help="number of time steps (default: %(default)s)")
parser.add_argument('--Nx', dest='Nxinner', action='store', default=2048, type=int, help="number of cells (default: %(default)s)")
parser.add_argument('--Nspecies', dest='Nspecies', action='store', default=1, type=int, help="number of species (default: %(default)s)")
parser.add_argument('--presolve', dest='presolve', action='store_true', help="Perform Poisson solve for E_x before timestepping (default: %(default)s)")
parser.add_argument('--theta', dest='theta', action='store', default=0.5, type=float, help="theta (default: %(default)s)")
parser.add_argument('--tol', dest='tol', action='store', default=1e-6, type=float, help="tol (default: %(default)s)")
parser.add_argument('--atol', dest='atol', action='store', default=1e-15, type=float, help="atol (default: %(default)s)")
parser.add_argument('--nppc', dest='nppc', action='store', default=64, type=int, help="number of particles per cell and species (default: %(default)s)")
parser.add_argument('--outputdir', dest='outputdir', action='store', default="run1", type=str, help="Output directory (default: %(default)s)")
parser.add_argument('--rescale_dx', dest='rescale_dx', action='store', default=1.0, type=float, help="Spatial resolution improvement (default: %(default)s)")
parser.add_argument('--rescale_dt', dest='rescale_dt', action='store', default=1.0, type=float, help="Temporal resolution improvement (default: %(default)s)")
parser.add_argument('--outputsteps', dest='outputsteps', action='store', default=1, type=int, help="Timesteps between field outputs (default: %(default)s)")
parser.add_argument('--particlesteps', dest='particlesteps', action='store', default=1000, type=int, help="Timesteps between particle outputs (default: %(default)s)")

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
            'build_A':0.,
            'build_b':0.,
            'current':0.,
            'density':0.,
            'energies':0.,
            'gmres':0.,
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
           and 3 ndarray(Np,float) for vx in cm/s

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

        if args.vbeam == 0.:
            self.vx = np.random.normal(scale=vthe, size=Np)
        else:
            self.vx = np.zeros(Np)
            Lx = args.Nxinner*dx
            for i in range(Np):
                self.vx[i] = (-1.)**i * args.vbeam*c + np.random.normal(scale=vthe, size=1) + args.deltav*c * np.sin(2.*np.pi*args.m*self.x[i]/Lx)

    def uncenter(self):
        """
        Pushes particle positions back from t=0 to t^-1/2

        Reads: x ndarray(Np,float): position in cm
               vx ndarray(Np,float): velocities in cm/s

        Changes: x to ndarray(Np,float) for oldx in cm

        Returns: nothing
        """
        Np = len(self.x)
        assert len(self.vx) == Np

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
               vx ndarray(Np,float): velocity in cm/s

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

    def newv(self, Elocalx):
        """
        Update velocity v from t^n to t^n+1 using E^n+theta at position x at t^n+1/2

        Elocalx ndarray(Np,float): electric field in statV/cm

        Reads: x ndarray(Np,float): position in cm
               vx ndarray(Np,float): velocity in cm/s

        Changes: ndarray(Np,float) for new positions as floats in cm, periodically
                 wrapped around if the particle leaves the (interior of) the domain

        Return: nothing
        """

        timers.tic("newv")

        Np = len(self.x)
        assert len(self.vx) == Np
        assert len(Elocalx) == Np

        beta = (self.q*dt)/(2.*self.m)

        vrx = self.vx + beta*Elocalx

        self.vx = 2.*vrx - self.vx

        timers.toc("newv")

    def get_localE(self, Ex):
        """
        Interpolate electric field from nodes to position of each particle

        Ex ndarray(Nxouter,float): Electric field in statV/cm

        Reads: x ndarray(Np,float): Particle position in cm

        Returns: 3 ndarray(Np,float) for the local electric field in statV/cm at
                 each particle location
        """
        timers.tic("localE")

        assert len(Ex) == args.Nxouter

        cellx = np.array((self.x/dx), dtype=int)

        wxfl = W_CIC(self.x, cellx-1, 0.0)
        wxfc = W_CIC(self.x, cellx  , 0.0)
        wxfr = W_CIC(self.x, cellx+1, 0.0)

        Elocalx = wxfl * Ex[cellx-1] + wxfc * Ex[cellx  ] + wxfr * Ex[cellx+1]

        timers.toc("localE")
        return Elocalx

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
def ecsim_current(Nx, x, vx, q,macro):
    """
    Current deposition of \hat{J}_{sg}

    Nx int: Nxouter
    x ndarray(Np,float): particle position in cm at t^n+1/2
    vx ndarray(Np,float): particle velocities in cm/s at t^n

    Returns: 3 ndarray(Nxouter,float) for currents at t^n in statC/cm^2/s already
             reduced to the interior of the periodic domain
    FIXME: is that current really at t^n+1/2?
    FIXME: I assume j is at integer options same as E?
    """
    jx = np.zeros(Nx)

    q = q * macro / dx**3

    Np = len(x)
    assert len(vx) == Np

    for p in range(Np):
        cellx = int(x[p]/dx)

        # for i in range(-2,3):
        for i in range(-1,2):
            w = W_CIC_scalar_numba(x[p], cellx+i, 0.0)
            jx[cellx+i] += q * w * vx[p]


    jx = charge_boundary(jx)
    return jx

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
def naive_numba_mass_matrix(Nx, x, q, m, macro, M00):
    """
    Compute mass matrices that contain the reaction of the magnetised plasma to
    the electric field

    Nx int: Nxouter
    x ndarray(Np,float): particle positions in cm at t^n+1/2

    Returns: 9 ndarray((Nxouter,Nxouter),float) for the mass matrices. The units of each
             entry are 1/seconds since we have to multiply an electric field in
             statV/cm and get a current in statC/cm^2/s
    """
    #assert len(x) == Np
    assert M00.shape[0] == Nx
    assert M00.shape[1] == Nx

    Vg = dx**3
    beta = (q*dt)/(2.*m*c) # dimension l^1/2 t m^-1/2

    for p in range(np.int64(len(x))):
        cellx = int(x[p]/dx)
        for i in range(-2,3):
            for j in range(-2,3):
                w = beta/Vg * q*macro*c * W_CIC_scalar_numba(x[p], cellx+i, 0.0) * W_CIC_scalar_numba(x[p], cellx+j, 0.0)
                M00[cellx+i,cellx+j] += w

    M00 = matrix_boundary(M00)

def init_fields():
    """
    Create fields at t=0

    Returns: ndarray(Nxouter,float) for Ex in statC/cm
    """

    Ex = np.zeros(args.Nxouter )
    return Ex

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

def maxwell(Ex, jhatx, M00, theta):
    """
    Solve linear problem to get updated electric field

    Ex ndarray(Nxouter,float): electric field at t^n in statV/cm

    jhatx ndarray(Nxouter,float): explicit current in statC/cm^2/s

    M00 ndarray((Nxouter,Nxouter),float): mass matricies in 1/s

    Returns: ndarray(Nxouter,float) for Ex at t^n+1 in statV/cm
    """

    timers.tic("maxwell")

    assert len(Ex) == args.Nxouter

    timers.tic("build_b")
    # our right hand side has units of Gauss/(cm/s * s) = (statV/cm)/cm = m^1/2 l^-3/2 t^-1
    b = np.zeros(args.Nxinner)

    jx = jhatx + (M00@Ex)*(1.-theta)

    for m in range(0*args.Nxinner,1*args.Nxinner):
        i = m - 0*args.Nxinner
        b[m] = (- Ex[i+Ng] / (c*dt) + 4*np.pi/c * jx[i+Ng])

    timers.toc("build_b")

    # are are constructing a matrix M for the equation
    # M.Ex = b where
    # Ex has units of      m^1/2 l^-1/2 t^-1
    # and b has units of  m^1/2 l^-3/2 t^-1
    # this implies that all entries in M should have units 1/l
    timers.tic("build_A")
    A = np.zeros( (args.Nxinner,args.Nxinner) )

    A[0*args.Nxinner:1*args.Nxinner, 0*args.Nxinner:1*args.Nxinner] -= (M00[Ng:-Ng, Ng:-Ng] * 4*np.pi/c * theta)
    for m in range(0*args.Nxinner,1*args.Nxinner):
        i = m - 0*args.Nxinner
        A[m,i] += (-1./(c*dt))

    timers.toc("build_A")

    timers.tic("gmres")
    info = 1
    tol = args.tol
    atol = args.atol
    while info > 0 and tol < 1e-2 and atol < 1e-2:
        Fnext,info = gmres(A,b, tol=tol/c, atol=atol, x0=Ex[Ng:-Ng])
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
    newEx[Ng:-Ng] = Fnext[0*args.Nxinner:1*args.Nxinner]

    newEx = field_boundary(newEx)

    timers.toc("maxwell")

    return newEx

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
    h5f.close()

def save_energies(energies, t, Ex, species):
    """
    Save enery at time t^n+1

    energies: open file pointer
    t: int time step we just completed
    Ex: electric field at t^n+1
    species: list of species that each has: position t^n+1/2 and velocity t^n+1
    """
    timers.tic("energies")

    eEx = (np.sum(Ex[Ng:-Ng]**2) * dx**3 / (8.*np.pi))
    eE  = eEx
    eFields = eE
    eParticles = np.zeros(len(species))
    for i,s in enumerate(species):
        E = 0.5 * s.m * s.macro * np.sum(s.vx**2)
        eParticles[i] = E

    energies.write(str((t*dt*args.wpe)))                   # 1
    energies.write(" "+str((np.sum(eParticles)+eFields)))  # 2
    energies.write(" "+str(np.sum(eParticles)))            # 3
    energies.write(" "+str(eFields))                       # 4
    energies.write(" "+str(0.))                            # 5
    energies.write(" "+str(eE))                            # 6
    energies.write(" "+str(0.))                            # 7
    energies.write(" "+str(0.))
    energies.write(" "+str(0.))                            # 9
    energies.write(" "+str(eEx))
    energies.write(" "+str(0.))                            # 11
    energies.write(" "+str(0.))                            # 12
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
sys.stderr.write("lD = "+str(lD)+"cm\n")
de = c/args.wpe
sys.stderr.write("de = "+str(de)+"cm\n")

#dx = .66713366*lD
dx = lD/args.rescale_dx
sys.stderr.write("dx = "+str(dx)+"cm\n")
sys.stderr.write("   = "+str(dx/lD)+"l_D\n")
sys.stderr.write("   = "+str(dx/de)+"d_e\n")


sys.stderr.write("Lx = "+str(dx*args.Nxinner)+"cm\n")
sys.stderr.write("   = "+str(dx*args.Nxinner/lD)+"l_D\n")
sys.stderr.write("   = "+str(dx*args.Nxinner/de)+"d_e\n")

dt_wpe = 1./args.wpe
dt_part = dx/(4.*vthe)
dt_yee = 0.9999 * dx/(np.sqrt(2.)*c)
dt = min(dt_wpe,dt_part,dt_yee) / args.rescale_dt

sys.stderr.write("dt = "+str(dt)+"s\n")
sys.stderr.write("wpe*dt = "+str((args.wpe*dt))+"\n")
sys.stderr.write("c*dt/dx = "+str((c*dt/dx))+"\n")

sys.stderr.write("T = "+str(dt*args.Nt)+"s\n")
sys.stderr.write("  = "+str(dt*args.Nt*args.wpe)+"wpe^-1\n")


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
    sys.stderr.write("Lx = "+str(dx*args.Nxinner/di)+"d_i\n")

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

t=0

Ex = init_fields()

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

# f.write("0 "+str(electrons.x[0])+" "+str(electrons.vx[0])+"\n")

# reset HDF5 output
h5f = h5py.File(args.outputdir+"/dings.h5", 'w')
cfg_grp = h5f.create_group("config")
cfg_grp.create_dataset("plasmafreq", data=np.array([args.wpe]))
cfg_grp.create_dataset("B0", data=np.array([0.,0.,0.]))
cfg_grp.create_dataset("units", data=np.array("cgs", dtype=h5py.string_dtype('ascii', 3)))
catfile = "lineout_raumschritte = 1\nmpzume="+str(args.mime)+"\nwidthbg="+str(vthe/c)+"\n"
utf8_type=h5py.string_dtype('utf-8', len(catfile))
catbytes = np.array(catfile.encode("utf-8"), dtype=utf8_type)
cfg_grp.create_dataset("cat-file", data=catbytes)
h5f.close()

saveh5(rho,  "/Timestep_0/rho/rhoL/rhoL[0]")
saveh5(rhoE, "/Timestep_0/rho/rhoE/rhoE")
saveh5(Ex, "/Timestep_0/felder/E/E[0]")
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
energies.flush()
save_energies(energies, 0, Ex, species)

energycheck = open(args.outputdir+"/energycheck.dat", 'w')
charge = open(args.outputdir+"/charge.dat", 'w')

momentum = open(args.outputdir+"/momentum.dat", "w")
p0 = 0.
px = np.zeros(args.Nt+1)
for s in species:
    p0 += s.m * s.macro * args.vbeam*c * len(s.vx)
    px[0] += s.m * s.macro * np.sum(s.vx)
if args.vbeam == 0.:
    for s in species:
        p0 += me * s.macro * args.vthe*c * len(s.vx)
momentum.write(str(0.)+" "+str(p0)+" "+str(px[0])+"\n")

timers.toc("init")
timers.tic("timestepping")

for t in range(1,args.Nt+1):
    timers.tic("misc1")
    # net current from all species
    jx = np.zeros(args.Nxouter)

    # total mass matrix for all species
    M00 = np.zeros( (args.Nxouter,args.Nxouter) )
    timers.toc("misc1")

    for s in species:
        # Update x from t^n-1/2 to t^n+1/2 using v at t^n
        s.newx()

        # compute current \hat{j}_sg using x at t^n+1/2 and B and v at t^n
        timers.tic('current')
        jhatx =  ecsim_current(args.Nxouter, s.x, s.vx, s.q,s.macro)
        timers.toc('current')

        # save for analysis
        # if args.outputsteps > 0 and t%args.outputsteps == 0:
        #     saveh5(jhatx, "/Timestep_"+str(t)+"/"+s.name+"/jhat[1]")

        # add to total current
        timers.tic("misc2")
        jx += jhatx
        timers.toc("misc2")

        # compute mass matrices
        timers.tic('mass_matrix')
        naive_numba_mass_matrix(args.Nxouter, s.x, s.q, s.m, s.macro, M00)
        timers.toc('mass_matrix')

    # save total current
    if args.outputsteps > 0 and t%args.outputsteps == 0:
        saveh5(jx,  "/Timestep_"+str(t)+"/rho/rhoL/rhoL[1]")

    # update E from t^n to t^n+1 using j at t^n+1/2 and the implicit theta algorithm
    newEx= maxwell(Ex, jx, M00, args.theta)

    # get electric field at t^n+1/2
    timers.tic("misc4")
    midEx = args.theta*newEx + (1.-args.theta)*Ex
    timers.toc("misc4")

    # save new fields
    if args.outputsteps > 0 and t%args.outputsteps == 0:
        saveh5(newEx, "/Timestep_"+str(t)+"/felder/E/E[0]")

    jxgbar = np.zeros(args.Nxouter)
    Ekinprev = 0.
    Ekinnext = 0.

    for s in species:
        # local electric field
        Elocalx = s.get_localE(midEx)

        # keep old velocity
        s.oldvx = s.vx.copy()

        # update particle velocities from t^n to t^n+1 using E at t^n+1
        s.newv(Elocalx)
        # check energy conservation
        s.vxbar = 0.5*(s.oldvx + s.vx)
        jxsgbar = ecsim_current(args.Nxouter, s.x, s.vxbar, s.q,s.macro)
        jxgbar += jxsgbar
        Ekinprev += 0.5*s.m*s.macro*np.sum(s.oldvx**2)
        Ekinnext += 0.5*s.m*s.macro*np.sum(s.vx**2)

    deltaEparticle = Ekinnext-Ekinprev
    deltaEfieldsEx = (np.sum(newEx[Ng:-Ng]**2)/(8.*np.pi) - np.sum(Ex[Ng:-Ng]**2)/(8.*np.pi)) * dx**3
    deltaEfields = deltaEfieldsEx
    work = np.sum(jxgbar[Ng:-Ng]*midEx[Ng:-Ng]) * dt * dx**3
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

    # update electric field to be at t^n+1
    timers.tic("misc6")
    Ex = newEx
    timers.toc("misc6")

    if args.particlesteps > 0 and t%args.particlesteps == 0:
        for s in species:
            save_particles(s, t)

    # save energy at t^n+1
    save_energies(energies, t, Ex, species)

    dispEx[t-1,:] = Ex[Ng:-Ng] * np.sin(np.pi*t/args.Nt)**2

    # f.write(str(t)+" "+str(electrons.x[0])+" "+str(electrons.vx[0])+"\n")

    # plt.scatter(electrons.x/dx, electrons.vx/c)
    # plt.show()

    for s in species:
        px[t] += s.m * s.macro * np.sum(s.vx)
    momentum.write(str(t*dt*args.wpe)+" "+str(p0)+" "+str(px[t])+"\n")

disp_kw = np.fft.fft2(dispEx)
power = np.abs(np.fft.fftshift(disp_kw.T)[:,args.Nt//2:])
save2d(power, "disp_E0_x.dat", scale=-1.)

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
print("Initalization:  "+str(timers.timers["init"]))
print("Time Stepping:  "+str(timers.timers["timestepping"]))
print("Mass Matrix:    "+str(timers.timers["mass_matrix"]))
print("Maxwell cost:   "+str(timers.timers["maxwell"]))
print("  build_b cost: "+str(timers.timers["build_b"]))
print("  build_A cost: "+str(timers.timers["build_A"]))
print("  GMRES cost:   "+str(timers.timers["gmres"]))
print("Current cost:   "+str(timers.timers["current"]))
print("Density cost:   "+str(timers.timers["density"]))
print("Local E:        "+str(timers.timers["localE"]))
print("New X:          "+str(timers.timers["newx"]))
print("New V:          "+str(timers.timers["newv"]))
print("Energies:       "+str(timers.timers["energies"]))
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
