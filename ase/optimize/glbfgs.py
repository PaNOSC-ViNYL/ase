import numpy as npy
import copy 
from ase.optimize import Optimizer
from ase.neb import *
class GLBFGS(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, dR=None,
                 memory=None, alpha=None, min='line'):
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

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

    def sign(self,w):
        if(w<0.0): return -1.0
        return 1.0

    def step(self, f):
        atoms = self.atoms
       # self.ni = atoms.nimages-2
       # try: atoms.imax
       # except: atoms.imax=0
       # if(not self.ni):atoms.imax=1
        self.r = npy.zeros((atoms.nimages * atoms.natoms, self.dim), 'd')
        for i in range(0, atoms.nimages):
            self.r = atoms.get_positions()
        self.f = npy.zeros((atoms.nimages * atoms.natoms, self.dim), 'd')
        for i in range(0, atoms.nimages):
            self.f = atoms.get_forces()
        try: self.start
        except:self.start=0
        if(not self.start):
            self.start = 1
            self.a = npy.zeros(self.memory+1, 'd')
            self.ptmp = atoms
            self.maxstep = npy.sqrt(self.maxstep * atoms.nimages)
            self.lbfgsinit = 0
        try: self.lbfgsinit
        except:self.lbfgsinit=0
        if(not self.lbfgsinit):
            self.lbfgsinit = 1
            self.Ho = npy.ones((atoms.nimages * atoms.natoms, self.dim), 'd')
            if (not self.min=='line'):self.Ho = self.Ho * self.alpha
            self.ITR = 1
            self.s = [1.]
            self.y = [1.]
            self.rho = [1.]
        else:
            a1 = abs (npy.vdot(self.f, self.f_old))
            a2 = npy.vdot(self.f_old, self.f_old)
            if(self.min=='line'):
                if(a1<=0.5* a2 and a2!=0):
                    reset_flag = 0
                else:
                    reset_flag = 1
            else:
                reset_flag = 0
            if(reset_flag==0):
                ITR = self.ITR#correctly generated
                if(ITR > self.memory):
                    self.s.pop(1)
                    self.y.pop(1)
                    self.rho.pop(1)
                    ITR=self.memory
                self.s.append(self.r - self.r_old)#!!self.r is not updating
        # boundry cond
        #        for i in range(atoms.ni):
        #            if(method=='min'):i=0
        #            try:
        #                DBC(self.s[ITR][i],atoms.p[i].Box) #need to make matrix for box
        #            except: 
        #                print "Box not found."
        #             if(method=='min'):break
                self.y.append(-(self.f-self.f_old))
                self.rho.append(1/npy.vdot(self.y[ITR],self.s[ITR]))
                self.ITR += 1
            else:
#        print 'reset image',i
                self.ITR = 1
                self.s = [1.]
                self.y = [1.]
                self.rho = [1.]
       
        self.r_old = self.r.copy()
        self.f_old = self.f.copy()

        if(self.ITR <= self.memory):
            BOUND = self.ITR
        else:
            BOUND = self.memory
        q = -1.0*self.f
        for j in range(1,BOUND):
            k = (BOUND-j)
            self.a[k] = self.rho[k] * npy.vdot(self.s[k], q)
            q -= self.a[k] * self.y[k]
        d = self.Ho * q # no needed cause Ho is idenity matrix 
        for j in range(1,BOUND):
            B = self.rho[j] * npy.vdot(self.y[j], d)
            d= d + self.s[j] * (self.a[j] - B)

        d = -1.0 * d
        if(not self.min=='line'): self.d = d
        self.du = d / npy.sqrt(npy.vdot(d, d))
        if(self.min=='line'):

      # Finite difference step using temporary point
            self.ptmp.r = self.r.copy()
            self.ptmp.imax = atoms.imax
#   must be turned on if tanset is commented
#        self.ptmp.p[i].N=atoms.p[i].N.copy()

            self.ptmp.r += (self.du * self.dR)

      # put back in points
            for i in range(0,atoms.nimages):
                self.ptmp.images[i].r = self.ptmp.r[i-1]
                self.ptmp.images[i].f=self.ptmp.images[i].get_forces()

        else:
    # use the Hessian Matrix to predict the min
            if(abs(npy.sqrt(self.d * self.d).sum()) > self.maxstep):
                self.d = self.du * self.maxstep
            self.r += self.d
            atoms.set_positions(self.r)
        if(self.min=='line'):
     # force projections on a small stepa
#            if(method=='neb' or method=='str' ):
#                self.ptmp.tangentset()# if commented copy N above
#                self.ptmp.project()
      # put force in big vector
            self.ptmp.f = npy.zeros((atoms.nimages * atoms.natoms, self.dim), 'd')
            for i in range(0,atoms.nimages):
                self.ptmp.f[i-1] = self.ptmp.images[i].f

      # Decide how much to move along the line du
            Fp1=npy.vdot(self.f, self.du)
            Fp2=npy.vdot(self.ptmp.f, self.du)
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
            self.r += (self.du * RdR)

    # put back in points
        for i in range(0,atoms.nimages):
            atoms.images[i].r=self.r[i-1]
