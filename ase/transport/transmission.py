import numpy as npy
from numpy import linalg
from ase.transport.tools import dagger

class Transmission:
    def __init__(self, **kwargs):
        self.input_parameters = {'energies': None,
                                 'greenfunction': None,
                                 'selfenergies': None,
                                 'transmission': True,
                                 'eigenchannels': 0,
                                 'dos': True,
                                 'pdos': [],
                                 'verbose': False}

        self.set(**kwargs)
        self.initialized = False

    def set(self, **kwargs):
        p = self.input_parameters
        p.update(kwargs)
        self.uptodate = False
        self.initialized = False #XXX smarter logic needed!
        
    def initialize(self):
        p = self.input_parameters
        self.verbose = p['verbose']
        self.energies = p['energies']
        self.nepts = len(self.energies)
        self.selfenergies = p['selfenergies']
        self.greenfunction = p['greenfunction']
        if p['transmission']:
            self.T_e = npy.empty(self.nepts)
        if p['dos']:
            self.dos_e = npy.empty(self.nepts)
        if len(p['pdos']) != 0:
            self.pdos_ne = npy.empty((len(p['pdos']), self.nepts))
        if p['eigenchannels'] > 0:
            self.eigenchannels_ne = npy.empty((p['eigenchannels'], self.nepts))
       
        self.uptodate = False
        self.initialized = True
            
    def update(self):
        if not self.initialized:
            self.initialize()

        p = self.input_parameters
        for e in range(self.nepts):
            if self.verbose:
                print "%i out of %i" % (e + 1, self.nepts)
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
                
    def get_t_matrix(self, e):
        Ginv_mm = self.greenfunction(self.energies[e], inverse=True)
        lambda1_mm = self.selfenergies[0].get_lambda(self.energies[e])
        lambda2_mm = self.selfenergies[1].get_lambda(self.energies[e])
        a_mm = linalg.solve(Ginv_mm, lambda1_mm)
        b_mm = linalg.solve(dagger(Ginv_mm), lambda2_mm)
        T_mm = npy.dot(a_mm, b_mm)
        return T_mm 
        
    def calculate_transmission_and_eigenchannels(self, e):
        T_mm = self.get_t_matrix(e)
        nchan = len(self.eigenchannels_ne)
        #self.T_e[e] = npy.trace(T_mm).real
        t_n = linalg.eigvals(T_mm).real
        self.eigenchannels_ne[:, e] = npy.sort(t_n)[-nchan:]
        self.T_e[e] = npy.sum(t_n)
   
    def calculate_transmission(self,e):
        self.T_e[e] = npy.trace(self.get_t_matrix(e)).real

    def calculate_dos(self, e):
        self.dos_e[e] = self.greenfunction.dos(self.energies[e])
        
    def calculate_pdos(self, e):
        bfs = self.input_parameters['pdos']
        pdos = self.greenfunction.pdos(self.energies[e])
        self.pdos_ne[:, e] = npy.take(pdos, bfs)
       
    def get_left_channels(self, energy, n=1):
        g_s_ii = self.greenfunction(energy)
        lambda_l_ii = self.selfenergies[0].get_lambda(energy)
        lambda_r_ii = self.selfenergies[1].get_lambda(energy)

        if self.greenfunction.S is None:
            s_s_qsrt_ii = s_s_isqrt = npy.identity(len(g_s_ii))
        else:
            s_mm = self.greenfunction.S
            s_s_i, s_s_ii = linalg.eig(s_mm)
            s_s_i = npy.abs(s_s_i)
            s_s_sqrt_i = npy.sqrt(s_s_i)#sqrt of eigenvalues  
            s_s_sqrt_ii = npy.dot(s_s_ii * s_s_sqrt_i, dagger(s_s_ii))
            s_s_isqrt_ii = npy.dot(s_s_ii / s_s_sqrt_i, dagger(s_s_ii))

        lambdab_r_ii = npy.dot(npy.dot(s_s_isqrt_ii, lambda_r_ii),s_s_isqrt_ii)
        a_l_ii = npy.dot(npy.dot(g_s_ii, lambda_l_ii), dagger(g_s_ii))
        ab_l_ii = npy.dot(npy.dot(s_s_sqrt_ii, a_l_ii), s_s_sqrt_ii)
        lambda_i, u_ii = linalg.eig(ab_l_ii)
        ut_ii = npy.sqrt(lambda_i / (2.0 * npy.pi)) * u_ii
        m_ii = 2 * npy.pi * npy.dot(npy.dot(dagger(ut_ii), lambdab_r_ii),ut_ii)
        T_i,c_in = linalg.eig(m_ii)
        T_i = npy.abs(T_i)
        channels = npy.argsort(-T_i)[:n]
        c_in = npy.take(c_in, channels, axis=1)
        T_n = npy.take(T_i, channels)
        v_in = npy.dot(npy.dot(s_s_isqrt_ii, ut_ii), c_in)

        return T_n, v_in
