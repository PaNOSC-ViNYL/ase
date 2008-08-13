import numpy as npy
import numpy.linalg as linalg
from ase.transport.tools import lambda_from_self_energy, dagger, tri2full

class GreensFunction:
    
    def __init__(self, **kwargs):
        
        self.input_parameters = {'energy' : None,
                                 'h_mm' : None,
                                 's_mm' : None,
                                 'selfenergies' : [],
                                 'eta' : 1.0e-4}

        self.initialized = False
        self.uptodate = False
        self.set(**kwargs)

    def initialize(self):
        p = self.input_parameters
        self.h_mm = p['h_mm']
        self.s_mm = p['s_mm']
        self.nbf = len(self.h_mm)
        self.energy = p['energy']
        self.selfenergies = p['selfenergies']
        self.eta = p['eta']
        #inverse of the scattering green's function
        self.gf_inv_mm = npy.empty((self.nbf,self.nbf),npy.complex) 
        self.initialized = True
        
    def set(self, **kwargs):
        p = self.input_parameters
        p.update(kwargs)
        self.energy = p['energy']
        self.uptodate = False

    def set_energy(self,energy):
        if self.energy == None:
            self.energy = energy
        elif abs(energy-self.energy) > 1.0e-12:
            self.energy = energy
            self.uptodate = False
        
    def update(self):
        """force and update""" 
        z =  self.energy + self.eta * 1.0j
        #update selfenergies
        sigmas_mm = []
        for selfenergy in self.selfenergies:
            selfenergy.set_energy(self.energy)
            if not selfenergy.uptodate:
                selfenergy.update()
            sigmas_mm.append(selfenergy.sigma_mm)
        
        #self.gf_inv_mm[:] = z*self.s_mm - self.h_mm - npy.sum(sigmas_mm,axis=0)
        self.gf_inv_mm[:] = self.s_mm
        self.gf_inv_mm *= z
        self.gf_inv_mm -= self.h_mm
        self.gf_inv_mm -= npy.sum(sigmas_mm,axis=0)
        self.uptodate = True

    def get_inv_matrix(self):
        if not self.uptodate:
            self.update()
        return self.gf_inv_mm

    def get_matrix(self):
        if not self.uptodate:
            self.update()
        return linalg.inv(self.g_inv_mm)

