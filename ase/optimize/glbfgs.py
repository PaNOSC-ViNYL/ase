import numpy as npy
import copy 
from ase.optimize import Optimizer
from ase.neb import *
class GLBFGS(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', maxstep=None, dR=None,
                 memory=None, alpha=None, min='line'):
        Optimizer.__init__(self, atoms, restart, logfile)

        if maxstep is not None:
            self.maxstep = maxstep
        if dR is not None:
            self.dR = dR
        if memory is not None:
            self.memory = memory
        if alpha is not None:
            self.alpha = alpha
        if min is not None:
            self.min = min
        self.min = min

    def initialize(self):
        self.H = None
        self.maxstep = 0.2
        self.dR = 0.0001
        self.memory = 25
        self.alpha = 0.05
        self.dim = 3
#        self.min = 'line'

    def sign(self,a):
        if(a<0.0): return -1.0
        return 1.0

    def step(self, f):
        atoms = self.atoms
        self.ni = atoms.nimages-2
        try: atoms.imax
        except: atoms.imax=0
        if(not self.ni):atoms.imax=1

        atoms.r = npy.zeros((self.ni, atoms.natoms, self.dim), 'd')
        for i in range(1, atoms.nimages-1):
            atoms.r[i-1] = atoms.images[i].get_positions()
        atoms.f = npy.zeros((self.ni, atoms.natoms, self.dim), 'd')
        for i in range(1, atoms.nimages-1):
            atoms.f[i-1] = atoms.images[i].get_forces()
        try: atoms.start
        except:atoms.start=0
        if(not atoms.start):
            atoms.start = 1
            atoms.a = npy.zeros(self.memory+1, 'd')
#            self.ptmp = copy.deepcopy(atoms)
            self.ptmp = atoms
            self.maxstep = npy.sqrt(self.maxstep * self.ni)
            atoms.lbfgsinit = 0
        try: atoms.lbfgsinit
        except:atoms.lbfgsinit=0
        if(not atoms.lbfgsinit):
            atoms.lbfgsinit = 1
            atoms.Ho = npy.ones((self.ni, atoms.natoms, self.dim), 'd')
            atoms.ITR = 1
            atoms.s = [1.]
            atoms.y = [1.]
            atoms.rho = [1.]
        else:
            a1 = abs (npy.vdot(atoms.f, atoms.f_old))
            a2 = npy.vdot(atoms.f_old, atoms.f_old)
            if(self.min=='line'):
                if(a1<=0.5* a2 and a2!=0):
                    reset_flag = 0
                else:
                    reset_flag = 1
            else:
                reset_flag = 0
            if(reset_flag==0):
                ITR = atoms.ITR
                if(ITR > self.memory):
                    atoms.s.pop(1)
                    atoms.y.pop(1)
                    atoms.rho.pop(1)
                    ITR=self.memory
                atoms.s.append(atoms.r - atoms.r_old)
        # boundry cond
        #        for i in range(atoms.ni):
        #            if(method=='min'):i=0
        #            try:
        #                DBC(atoms.s[ITR][i],atoms.p[i].Box) #need to make matrix for box
        #            except: 
        #                print "Box not found."
        #             if(method=='min'):break
                atoms.y.append(-(atoms.f-atoms.f_old))
                atoms.rho.append(1/npy.vdot(atoms.y[ITR],atoms.s[ITR]))
                atoms.ITR += 1
            else:
#        print 'reset image',i
                atoms.ITR = 1
                atoms.s = [1.]
                atoms.y = [1.]
                atoms.rho = [1.]
        atoms.r_old = atoms.r.copy()
        atoms.f_old = atoms.f.copy()

        if(atoms.ITR <= self.memory):
            BOUND = atoms.ITR
        else:
            BOUND = self.memory
        q = -1.0*atoms.f
        for j in range(1,BOUND):
            k = (BOUND-j)
            atoms.a[k] = atoms.rho[k] * npy.vdot(atoms.s[k], q)
            q -= atoms.a[k] * atoms.y[k]

        d = atoms.Ho * q # no needed cause Ho is idenity matrix

        for j in range(1,BOUND):
            B = atoms.rho[j] * npy.vdot(atoms.y[j], d)
            d= d + atoms.s[j] * (atoms.a[j] - B)

        d = -1.0 * d
        if(not self.min=='line'): atoms.d = d
        atoms.du = d / npy.sqrt(npy.vdot(d, d))
        if(self.min=='line'):

      # Finite difference step using temporary point
            self.ptmp.r = atoms.r.copy()
            self.ptmp.imax = atoms.imax
#   must be turned on if tanset is commented
#        self.ptmp.p[i].N=atoms.p[i].N.copy()

            self.ptmp.r += (atoms.du * self.dR)

      # put back in points
            for i in range(1,atoms.nimages-1):
                self.ptmp.images[i].r = self.ptmp.r[i-1]
                self.ptmp.images[i].f=self.ptmp.images[i].get_forces()

        else:
    # use the Hessian Matrix to predict the min
            if(abs(npy.sqrt(atoms.d * atoms.d).sum()) > self.maxstep):
                atoms.d = atoms.du * self.maxstep
            atoms.r += atoms.d

        if(self.min=='line'):
     # force projections on a small stepa
#            if(method=='neb' or method=='str' ):
#                self.ptmp.tangentset()# if commented copy N above
#                self.ptmp.project()
      # put force in big vector
            self.ptmp.f = npy.zeros((self.ni, atoms.natoms, self.dim), 'd')
            for i in range(1,atoms.nimages-1):
                self.ptmp.f[i-1] = self.ptmp.images[i].f

      # Decide how much to move along the line du
            Fp1=npy.vdot(atoms.f, atoms.du)
            Fp2=npy.vdot(self.ptmp.f, atoms.du)
            CR=(Fp1 - Fp2) / self.dR
            if(CR < 0.0):
                print "negcurve"
                RdR = self.maxstep
#        RdR=Fp1*0.1
                if(abs(RdR) > self.maxstep):
                    RdR = self.sign(RdR) * self.maxstep
            else:
                Fp = (Fp1 + Fp2) * 0.5
                RdR = Fp / CR
                if(abs(RdR) > self.maxstep):
                    RdR = self.sign(RdR) * self.maxstep
#          print "max"
                else:
                    RdR += self.dR * 0.5
      # move to the next space
            atoms.r += (atoms.du * RdR)

    # put back in points
        for i in range(1,atoms.nimages-1):
            atoms.images[i].r=atoms.r[i-1]
