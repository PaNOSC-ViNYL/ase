import numpy as npy
from numpy import linalg
from ase.transport.tools import dagger

class Transmission:

    def __init__(self,**kwargs):
        
        self.input_parameters = {'energies' : None,
                                 'greensfunction' : None,
                                 'selfenergies' : None,
                                 'transmission' : True,
                                 'eigenchannels' : 0,
                                 'dos' : True,
                                 'pdos' : []}

        self.set(**kwargs)
        self.initialized = False

    def set(self,**kwargs):
        p = self.input_parameters
        p.update(kwargs)
        self.uptodate = False
        self.initialized = False #XXX smarter logic needed!
        
    def initialize(self):
        p = self.input_parameters
        self.energies = p['energies']
        self.nepts = len(self.energies)
        self.selfenergies = p['selfenergies']
        self.greensfunction = p['greensfunction']
        if p['transmission']:
            self.T_e = npy.empty(self.nepts)
        if p['dos']:
            self.dos_e = npy.empty(self.nepts)
        if len(p['pdos'])!=0:
            self.pdos_ne = npy.empty((len(p['pdos']),self.nepts))
        if p['eigenchannels']>0:
            self.eigenchannels_ne = npy.empty((p['eigenchannels'],self.nepts))
       
        self.uptodate = False
        self.initialized = True
            
    def update(self):
        if not self.initialized:
            self.initialize()

        p = self.input_parameters
        for e in range(self.nepts):
            if p['transmission']:
                if p['eigenchannels'] > 0:
                    self.calculate_transmission_and_eigenchannels(e)
                else:
                    self.calculate_transmission(e)
            if p['dos']:
                self.calculate_dos(e)
            if p['pdos'] != []:
                self.calculate_pdos(e)
        self.uptodate = True
                
    def get_t_matrix(self,e):
        gf = self.greensfunction
        energy = self.energies[e]
        gf.set_energy(energy)
        if not gf.uptodate:
            gf.update()
        lambda1_mm = self.selfenergies[0].get_lambda_matrix()
        lambda2_mm = self.selfenergies[1].get_lambda_matrix()
        a_mm = linalg.solve(gf.gf_inv_mm,lambda1_mm)
        b_mm = linalg.solve(dagger(gf.gf_inv_mm),lambda2_mm)
        T_mm = npy.dot(a_mm,b_mm)
        return T_mm 
        
    def calculate_transmission_and_eigenchannels(self,e):
        T_mm = self.get_t_matrix(e)
        nchan = len(self.eigenchannels_ne)
        #self.T_e[e] = npy.trace(T_mm).real
        t_n = linalg.eigvals(T_mm).real
        self.eigenchannels_ne[:,e] = npy.sort(t_n)[-nchan:]
        self.T_e[e] = npy.sum(t_n)
   
    def calculate_transmission(self,e):
        self.T_e[e] = npy.trace(self.get_t_matrix(e)).real

    def calculate_dos(self,e):
        gf = self.greensfunction
        energy = self.energies[e]
        gf.set_energy(energy)
        if not gf.uptodate:
            gf.update()
        gfs_mm = linalg.solve(gf.gf_inv_mm, gf.s_mm)
        self.dos_e[e] = -1.0 / npy.pi * npy.trace(gfs_mm.imag)
        
    def calculate_pdos(self,e):
        bfs = self.input_parameters['pdos']
        gf = self.greensfunction
        s_mm = gf.s_mm
        energy = self.energies[e]
        gf.set_energy(energy)
        if not gf.uptodate:
            gf.update()
        sgfs_mm = npy.dot(s_mm,linalg.solve(gf.gf_inv_mm,s_mm))
        a_m = (npy.diagonal(sgfs_mm) * (1.0 / npy.diagonal(s_mm))).imag
        a_m *= -1.0 / npy.pi
        self.pdos_ne[:,e] = npy.take(a_m,bfs) 
        
        
#    def get_transmission(self):
#        self.input_parameters['transmission'] = True
#        if not self.uptodate:
#            self.update()
#            self.uptodate = True
#        return self.T_e

#    def get_dos(self):
#        self.input_parameters['dos'] = True
#        if not self.uptodate:
#            self.update()
#            self.uptodate = True
#        return self.dos_e
   
#    def get_eigenchannels(self):
#        if not self.uptodate:
#            self.update()
#            self.uptodate = True
#        return self.eigenchannels_ne
