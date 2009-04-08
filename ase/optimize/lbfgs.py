import numpy as npy
import copy 
from ase.optimize import Optimizer
from ase.neb import *
import random
#out = open('forces_bfgs','w')
#out2 = open('step_bfgs','w')
#out3 = open('displacement_bfgs','w')
#out4 = open('dpm_tmp_bfgs','w')
class LBFGS(Optimizer):
    def __init__(self, atoms, restart='restart_info', logfile='-', trajectory=None,
                 maxstep=0.2, dR=0.1,
                 memory=25, alpha=0.05, method='line'):
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        if maxstep > 1.0:
            raise ValueError(
                             'Wanna fly? I know the calculation is too slow. ' +
                             'But you have to follow the rules.\n'+
                             '            The maximum step size %.1f' % maxstep +' is too big! \n'+
                             '            Try to set the maximum step size below 0.2.')
        self.maxstep = maxstep
        self.dR = dR
        self.memory = memory
        self.alpha = alpha
        self.method = method

    def initialize(self):
        self.lbfgsinit = 0
        self.ITR = 1
        self.s = [1.]
        self.y = [1.]
        self.rho = [1.]

    def sign(self,w):
        if(w<0.0): return -1.0
        return 1.0

    def read(self):
        self.lbfgsinit, self.ITR, self.s, self.y, self.rho, self.r_old, self.f_old = self.load()

    def step(self, f):
        #print f
        for i in range(len(f)):
            if(npy.vdot(f[i],f[i])>0):
                f[i][0]+=random.uniform(10e-4,10e-3)*(-1)**random.randint(0,1)
                f[i][1]+=random.uniform(10e-4,10e-3)*(-1)**random.randint(0,1)
                f[i][2]+=random.uniform(10e-4,10e-3)*(-1)**random.randint(0,1)
        #print f
        self.r = self.atoms.get_positions()
        self.update(self.r, f)

        if(self.ITR <= self.memory):
            BOUND = self.ITR
        else:
            BOUND = self.memory
        q = -1.0*f
        for j in range(1,BOUND):
            k = (BOUND-j)
            self.a[k] = self.rho[k] * npy.vdot(self.s[k], q)
            q -= self.a[k] * self.y[k]
        d = self.Ho * q # no needed cause Ho is idenity matrix 
        for j in range(1,BOUND):
            B = self.rho[j] * npy.vdot(self.y[j], d)
            d= d + self.s[j] * (self.a[j] - B)
        d = -1.0 * d

        if(not self.method=='line'): self.d = d
        self.du = d / npy.sqrt(npy.vdot(d, d))

        if(self.method=='line'):
        # Finite difference step using temporary point
            tmp_r = self.r.copy()
            tmp_r += (self.du * self.dR)
            self.tmp.set_positions(tmp_r)
        # Decide how much to move along the line du
            Fp1=npy.vdot(f, self.du)
            Fp2=npy.vdot(self.tmp.get_forces(), self.du)
            CR=(Fp1 - Fp2) / self.dR
            #RdR = Fp1*0.1
            if(CR < 0.0):
                print "negcurve"
                RdR = self.maxstep
                if(abs(RdR) > self.maxstep):
                    RdR = self.sign(RdR) * self.maxstep
            else:
                Fp = (Fp1 + Fp2) * 0.5
                #Fp = (Fp1 + Fp2) * 0.3
                RdR = Fp / CR 
                if(abs(RdR) > self.maxstep):
                    RdR = self.sign(RdR) * self.maxstep
                else:
                    RdR += self.dR * 0.5
                    #RdR += self.dR * 0.3
            self.d = self.du * RdR####
        else:
        # use the Hessian Matrix to predict the min
            if(abs(npy.sqrt(npy.vdot(self.d, self.d).sum())) > self.maxstep):
                self.d = self.du * self.maxstep####
        self.r_old = self.r.copy()
        self.f_old = f.copy()
        self.r += self.d####
        self.atoms.set_positions(self.r)

    def update(self, r, f):
        self.start = 1
        self.a = npy.zeros(self.memory+1, 'd')
        self.tmp = self.atoms
        self.Ho = npy.ones((npy.shape(r)[0], 3), 'd')
        if (not self.method=='line'):self.Ho = self.Ho * self.alpha
        if(not self.lbfgsinit):
            self.lbfgsinit = 1
        else:
            a1 = abs (npy.vdot(f, self.f_old))
            a2 = npy.vdot(self.f_old, self.f_old)
            if(self.method=='line'):
                if(a1<=0.5* a2 and a2!=0):
                    reset_flag = 0
                else:
                    reset_flag = 1
            else:
                reset_flag = 0
            if(reset_flag==0):
                ITR = self.ITR
                if(ITR > self.memory):
                    self.s.pop(1)
                    self.y.pop(1)
                    self.rho.pop(1)
                    ITR=self.memory
                self.s.append(r - self.r_old)
                self.y.append(-(f-self.f_old))
                self.rho.append(1/npy.vdot(self.y[ITR],self.s[ITR]))
                self.ITR += 1
            else:
                self.ITR = 1
                self.s = [1.]
                self.y = [1.]
                self.rho = [1.]
            self.dump((self.lbfgsinit, self.ITR, self.s, self.y, self.rho, self.r_old, self.f_old))
