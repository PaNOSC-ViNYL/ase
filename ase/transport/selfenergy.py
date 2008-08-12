import numpy as npy
import numpy.linalg as linalg
from ase.transport.tools import dagger, tri2full, lambda_from_self_energy

class LeadSelfEnergy:
    
    def __init__(self,hs_dii,hs_dij,hs_dim,eta=1.0e-4,energy=None):
        self.h_ii = hs_dii[0] #onsite principal layer
        self.s_ii = hs_dii[1] 
        self.h_ij = hs_dij[0] #coupling between principal layers
        self.s_ij = hs_dij[1]
        self.h_im = hs_dim[0] #coupling to the central region
        self.s_im = hs_dim[1]
        self.eta = eta
        self.nbf = self.h_im.shape[1]
        self.sigma_mm = npy.empty((self.nbf,self.nbf), npy.complex)
        self.conv = 1.0e-3
        self.energy = None

    def get_hs_2pl_matrix(self):
        nbf = self.nbf
        h2pl = npy.zeros((2 * nbf,2 * nbf),npy.complex)
        s2pl = npy.zeros((2 * nbf,2 * nbf),npy.complex)
        h2pl[:nbf,:nbf] = self.h_ii
        h2pl[nbf:,nbf:] = self.h_ii
        h2pl[:nbf,nbf:2*nbf] = self.h_ij
        tri2full(h2pl,'U')

        s2pl[:nbf,:nbf] = self.s_ii
        s2pl[nbf:,nbf:] = self.s_ii
        s2pl[:nbf,nbf:2*nbf] = self.s_ij
        tri2full(s2pl,'U')

        return h2pl,s2pl

    def set_energy(self,energy):
        if self.energy == None:
            self.energy = energy
            self.uptodate = False
        elif abs(self.energy-energy) > 1.0e-12:
            self.energy = energy
            self.uptodate = False
    
    def get_ginv(self):
        """The inverse of the retarded surface GF""" 
        h_ii = self.h_ii
        s_ii = self.s_ii
        h_ij = self.h_ij
        s_ij = self.s_ij
        z = self.energy + self.eta * 1.0j
        
        v_00 = z * dagger(s_ii) - dagger(h_ii)
        v_11 = v_00.copy()
        v_10 = z * s_ij - h_ij
        v_01 = z * dagger(s_ij) - dagger(h_ij)

        delta = self.conv + 1
        n = 0
        while (delta > self.conv):
            a = linalg.solve(v_11,v_01)
            b = linalg.solve(v_11,v_10)
            v_01_dot_b = npy.dot(v_01,b)
            v_00 -= v_01_dot_b
            v_11 -= npy.dot(v_10,a) 
            v_11 -= v_01_dot_b
            v_01 = -npy.dot(v_01,a)
            v_10 = -npy.dot(v_10,b)
        
            delta = npy.abs(v_01).max()
            n += 1

        return v_00

    def update(self):
        """Force and update"""

        h_im = self.h_im
        s_im = self.s_im
        z = self.energy + self.eta * 1.j
            
        tau_im = z * s_im - h_im
        a_im = linalg.solve(self.get_ginv(),tau_im)
        tau_mi = z * dagger(s_im) - dagger(h_im)
        self.sigma_mm[:] = npy.dot(tau_mi,a_im)
        self.uptodate = True

    def get_matrix(self):
        if not self.uptodate:
            self.update()
        
        return self.sigma_mm
       
    def get_lambda_matrix(self):
        if not self.update:
            self.update()
        return lambda_from_self_energy(self.sigma_mm)
        
        
        
