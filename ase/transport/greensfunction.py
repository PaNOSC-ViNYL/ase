import numpy as npy

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
        return npy.linalg.inv(self.g_inv_mm)

## class GreenFunction:
##     """Equilibrium retarded Green function."""
    
##     def __init__(self, H, S=1, selfenergies=[], eta=1e-4):
##         self.H = H
##         self.S = S
##         self.selfenergies = selfenergies
##         self.eta = eta

##         setlf.energy = None
##         self.Ginv = npy.empty_like(H)

##     def __call__(self, energy, inverse=False):
##         if energy != self.energy:
##             self.energy = energy
##             self.Ginv[:] = (energy + self.eta * 1.j) * self.S - self.H

##             for selfenergy in self.selfenergies:
##                 selfenergy.set_energy(energy)
##                 if not selfenergy.uptodate:
##                     selfenergy.update()
##                 self.Ginv -= selfenergy.sigma_mm

##         if inverse:
##             return self.Ginv
##         else:
##             return npy.linalg.inv(self.Ginv)
